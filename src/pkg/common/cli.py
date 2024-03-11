import logging
import os
import signal


LOG_LEVEL = logging.INFO
LOG_FORMAT = "[%(asctime)s] %(message)s"
LOG_DATE_FORMAT = "%F %T"


def die(sig: int):
    signal.signal(sig, signal.SIG_DFL)
    os.kill(os.getpid(), sig)


def invoke_main(main, kwargs: dict, log: logging.Logger):
    """
    Starts a CLI main with given arguments, with standard signal-handling
    guards so the program behaves correctly in shell script.
    """
    try:
        logging.basicConfig(
            level=LOG_LEVEL, format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT
        )
        main(**kwargs)
    except KeyboardInterrupt:
        log.info("Canceled")
        die(signal.SIGINT)
    except BrokenPipeError:
        die(signal.SIGPIPE)
    except Exception:
        log.exception("Uncaught exception")
