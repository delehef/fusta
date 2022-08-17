#![allow(clippy::redundant_field_names)]
#[macro_use]
extern crate lazy_static;

use anyhow::{bail, Context, Result};
use clap::*;
use daemonize::*;
use log::*;
#[cfg(feature = "notifications")]
use notify_rust::Notification;
use simplelog::*;

pub mod fs;
use fs::*;

#[derive(Debug, Clone)]
struct RunEnvironment {
    mountpoint: std::path::PathBuf,
    created_mountpoint: bool,
}
fn main() -> Result<()> {
    let args =
        App::new("fusta")
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::UnifiedHelpMessage)
        .version(crate_version!())
        .author(crate_authors!())
        .arg(Arg::with_name("FASTA")
             .help("A (multi)FASTA file containing the sequences to mount")
             .required(true)
             .index(1))
        .arg(Arg::with_name("verbose")
             .short('v')
             .action(ArgAction::Count)
             .help("Sets the level of verbosity"))

        .arg(Arg::with_name("mountpoint")
             .short('o')
             .long("mountpoint")
             .help("Specifies the directory to use as mountpoint; it will be created if it does not exist")
             .takes_value(true))
        .arg(Arg::with_name("nodaemon")
             .short('D')
             .long("no-daemon")
             .help("Do not daemonize"))

    // Technical options
        .arg(Arg::with_name("max-cache")
             .short('C')
             .long("max-cache")
             .help("Set the maximum amount of memory to use to cache writes (MB)")
             .default_value("500")
             .takes_value(true))
        .arg(Arg::with_name("cache")
             .long("cache")
             .help("Use either mmap, fseek(2) or memory-backed cache to extract sequences from FASTA files. WARNING: memory caching use as much RAM as the size of the FASTA file should be available.")
             .possible_values(&["file", "mmap", "memory"])
             .default_value("mmap"))

    // Other options
        .arg(Arg::with_name("csv-separator")
             .short('S')
             .long("sep")
             .help("Set the separator to use in CSV files")
             .default_value(",")
             .takes_value(true))
        .get_matches();

    let log_level = match args.get_one::<u8>("verbose").copied().unwrap_or_default() {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        2 => LevelFilter::Trace,
        _ => LevelFilter::Trace,
    };
    let log_config = ConfigBuilder::new().build();
    let mut loggers: Vec<Box<dyn SharedLogger>> = vec![TermLogger::new(
        log_level,
        log_config.clone(),
        TerminalMode::Mixed,
        simplelog::ColorChoice::Auto,
    )];
    if !args.is_present("nodaemon") {
        let log_file_path = tempfile::Builder::new()
            .prefix("fusta-")
            .suffix(".log")
            .tempfile()
            .context("Unable to create a temporary file")?;

        println!(
            "Logs ({:?}) available in {}",
            log_level,
            log_file_path.path().display()
        );
        loggers.push(WriteLogger::new(log_level, log_config, log_file_path));
    }
    CombinedLogger::init(loggers).context("Unable to init logger")?;

    let fasta_file = value_t!(args, "FASTA", String)?;
    let mountpoint = value_t!(args, "mountpoint", String).unwrap_or(format!(
        "fusta-{}",
        std::path::Path::new(&fasta_file)
            .file_stem()
            .and_then(|s| s.to_str())
            .context(format!("{:?} is not a valid path", &fasta_file))?
    ));
    let fuse_options: Vec<fuser::MountOption> = vec![
        fuser::MountOption::FSName("FUSTA".to_string()),
        // fuser::MountOption::AutoUnmount,
        fuser::MountOption::DefaultPermissions,
    ];
    let settings = FustaSettings {
        cache: match args.value_of("cache").unwrap() {
            "mmap" => fs::Cache::Mmap,
            "file" => fs::Cache::File,
            "memory" => fs::Cache::RAM,
            _ => unreachable!(),
        },
        concretize_threshold: value_t!(args, "max-cache", usize).unwrap() * 1024 * 1024,
        csv_separator: value_t!(args, "csv-separator", String).unwrap(),
    };
    info!("Caching method:  {:#?}", settings.cache);

    let fs = FustaFS::new(settings, &fasta_file)?;

    let mut env = RunEnvironment {
        mountpoint: std::path::PathBuf::from(mountpoint),
        created_mountpoint: false,
    };
    if !env.mountpoint.exists() {
        std::fs::create_dir(&env.mountpoint)?;
        env.created_mountpoint = true;
    }
    if !env.mountpoint.is_dir() {
        bail!("mount point `{:?}` is not a directory", env.mountpoint);
    }
    if std::fs::read_dir(&env.mountpoint)?.take(1).count() != 0 {
        bail!("mount point {:?} is not empty.", env.mountpoint);
    }

    let umount_msg = if cfg!(target_os = "freebsd") || cfg!(target_os = "macos") {
        format!(
            "Please use `umount {:?}` to exit.",
            &env.mountpoint.canonicalize().unwrap()
        )
    } else {
        format!(
            "Please use `fusermount -u {0:?}` or `umount {0:?}` to exit.",
            &env.mountpoint.canonicalize().unwrap()
        )
    };

    info!("{}", &umount_msg);
    {
        ctrlc::set_handler(move || {
            error!("{}", umount_msg);
        })?;
    }

    if !args.is_present("nodaemon") {
        let pid_file = tempfile::Builder::new()
            .prefix("fusta-")
            .suffix(".pid")
            .tempfile()
            .context("Unable to create a temporary PID file")?;

        Daemonize::new()
            .pid_file(pid_file)
            .working_directory(std::env::current_dir().context("Unable to read current directory")?)
            .start()?;
    }

    #[cfg(feature = "notifications")]
    Notification::new()
        .summary("FUSTA")
        .body(&format!(
            "{} is now available in {:#?}",
            &fasta_file, &env.mountpoint
        ))
        .show()
        .unwrap();

    match fuser::mount2(fs, &env.mountpoint, &fuse_options) {
        Ok(()) => {}
        _ => {
            error!("Unable to mount the FUSE filesystem");
            std::process::exit(1);
        }
    }
    cleanup(&env)?;

    Ok(())
}

fn cleanup(env: &RunEnvironment) -> Result<()> {
    #[cfg(feature = "notifications")]
    Notification::new()
        .summary("FUSTA")
        .body("Successfully unmounted")
        .show()
        .unwrap();

    if env.created_mountpoint {
        #[cfg(feature = "notifications")]
        Notification::new()
            .summary("FUSTA")
            .body(&format!(
                "You can now safely remove the {:?} directory",
                env.mountpoint
            ))
            .show()
            .unwrap();
        warn!(
            "You can now safely remove the {:?} directory",
            env.mountpoint
        );
    }
    Ok(())
}
