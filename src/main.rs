#![allow(clippy::redundant_field_names)]
use std::sync::mpsc::channel;
use std::ffi::OsStr;

use daemonize::*;
use anyhow::{Result, bail, Context};
use log::*;
use simplelog::*;
use clap::*;

pub mod fs;
use fs::*;





#[derive(Debug, Clone)]
struct RunEnvironment {
    mountpoint: std::path::PathBuf,
    created_mountpoint: bool,
}
fn main() -> Result<()> {
    let args = App::new("fusta")
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::UnifiedHelpMessage)
        .version(crate_version!())
        .author(crate_authors!())
        .arg(Arg::with_name("FASTA")
             .help("A (multi)FASTA file containing the sequences to mount")
             .required(true)
             .index(1))
        .arg(Arg::with_name("mountpoint")
             .short("o")
             .long("mountpoint")
             .help("Specifies the directory to use as mountpoint")
             .default_value("fusta")
             .takes_value(true)
             .required(true))

        .arg(Arg::with_name("v")
             .short("v")
             .multiple(true)
             .help("Sets the level of verbosity"))
        .arg(Arg::with_name("daemon")
             .short("D")
             .long("daemon")
             .help("Launch in the background; will automatically quit when unmounted"))

        .arg(Arg::with_name("mmap")
             .short("M")
             .long("mmap")
             .help("Use mmap instead of seek to browse FASTA fragments. Faster, but memory hungry"))

        .arg(Arg::with_name("nonempty")
             .short("E")
             .long("non-empty")
             .help("Perform the mount even if the destination folder is not empty"))
        .get_matches();

    let log_level = match args.occurrences_of("v") {
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };
    let log_config = ConfigBuilder::new()
        .set_time_format_str("")
        .build();

    if args.is_present("daemon") {
        let mut log_file = std::env::temp_dir();
        log_file.push("fusta.log");

        WriteLogger::init(
            log_level,
            log_config,
            std::fs::File::create(&log_file).context("Unable to create log file")?
        ).context("Unable to initialize logger")?;

        let mut pid_file = std::env::temp_dir();
        pid_file.push("fusta.pid");

        println!("Logs {:?} available in {:?}", log_level, log_file);

        Daemonize::new()
            .pid_file(pid_file)
            .working_directory(std::env::current_dir().context("Unable to read current directory")?)
            .start()?;
    } else {
        TermLogger::init(log_level, log_config, TerminalMode::Mixed).context("Unable to initialize logger")?;
    }



    let fasta_file = value_t!(args, "FASTA", String)?;
    let mountpoint = value_t!(args, "mountpoint", String)?;
    let mut fuse_options: Vec<&OsStr> = vec![
        &OsStr::new("-o"), &OsStr::new("auto_unmount"),
        &OsStr::new("-o"), &OsStr::new("default_permissions"),
    ];
    let settings = FustaSettings {
        mmap: args.is_present("mmap"),
    };

    info!("Using MMAP:      {}", settings.mmap);

    let fs = FustaFS::new(settings, &fasta_file);
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

    if args.is_present("nonempty") {
        fuse_options.push(&OsStr::new("-o"));
        fuse_options.push(&OsStr::new("nonempty"));
    }
    if !args.is_present("nonempty") && std::fs::read_dir(&env.mountpoint)?.take(1).count() != 0 {
        bail!(
            "mount point {:?} is not empty, use the -E flag to mount in a non-empty directory.",
            env.mountpoint
        );
    }

    let (tx, rx) = channel();


    {
        let tx_ctrlc = tx.clone();
        ctrlc::set_handler(move || {
            info!("Ctrl-C received, exiting.");
            let _ = tx_ctrlc.send(());
        })?;
    }

    {
        let tx_fuse = tx.clone();
        let env = env.clone();
        let _ = std::thread::spawn(move || {
            match fuse::mount(fs, &env.mountpoint, &fuse_options) {
                Ok(()) => {},
                _      => { error!("Unable to mount the FUSE filesystem"); std::process::exit(1); }
            }
            let _ = tx_fuse.send(());
        });
    }

    rx.recv()?;
    cleanup(&env)?;

    Ok(())
}

fn cleanup(env: &RunEnvironment) -> Result<()> {
    if env.created_mountpoint {
        warn!("You can now safely remove the {:?} directory", env.mountpoint);
    }
    Ok(())
}
