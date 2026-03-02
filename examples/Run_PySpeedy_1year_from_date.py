import argparse
from datetime import datetime
from pathlib import Path

from pyspeedy import Speedy
from pyspeedy.callbacks import (
    DailyDiagnosticsExporter,
    DiagnosticCheck,
    RuntimeSummary,
    XarrayExporter,
)
from pyspeedy.config import load_config


def _add_calendar_year(value):
    try:
        return value.replace(year=value.year + 1)
    except ValueError:
        # Handle leap-day starts by falling back to Feb 28 the following year.
        return value.replace(year=value.year + 1, month=2, day=28)


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Run a one-year pySPEEDY integration from a user-specified start date."
    )
    parser.add_argument(
        "run_name",
        help="Name used to create a dedicated output subdirectory for this run.",
    )
    parser.add_argument(
        "start_date",
        help="Run start date in ISO format, for example 1980-01-01.",
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help="Optional parent output directory. Defaults to the configured run.output_dir.",
    )
    parser.add_argument(
        "--gwd",
        choices=("config", "on", "off"),
        default="config",
        help="Override the orographic gravity-wave-drag setting. Default: use config.",
    )
    return parser.parse_args()


args = _parse_args()
config = load_config()
run_config = config.run
model_config = config.model

start_date = datetime.fromisoformat(args.start_date)
end_date = _add_calendar_year(start_date)
summary_interval = 30 * model_config.nsteps
output_root = Path(args.output_root) if args.output_root is not None else Path(run_config.output_dir)
output_dir = output_root / args.run_name
output_dir.mkdir(parents=True, exist_ok=True)

model = Speedy(
    start_date=start_date,
    end_date=end_date,
)
if args.gwd == "on":
    model["orographic_gwd_enabled"] = True
elif args.gwd == "off":
    model["orographic_gwd_enabled"] = False
model.set_bc()

callbacks = [
    DiagnosticCheck(interval=run_config.diag_interval),
    RuntimeSummary(
        interval=summary_interval,
        verbose=run_config.verbose_output,
    ),
    XarrayExporter(
        interval=run_config.history_interval,
        output_dir=str(output_dir),
        variables=run_config.output_vars,
        verbose=run_config.verbose_output,
    ),
    DailyDiagnosticsExporter(
        output_dir=str(output_dir),
        verbose=run_config.verbose_output,
    ),
]
history_exporter = callbacks[2]

print(f"Run name: {args.run_name}")
print(f"Running from {start_date:%Y-%m-%d} to {end_date:%Y-%m-%d} (end date exclusive).")
print(f"Output directory: {output_dir}")
print(f"GWD enabled: {bool(model['orographic_gwd_enabled'])}")
print(f"Runtime summaries every {summary_interval // model_config.nsteps} days.")
print(history_exporter.variables)

model.run(callbacks=callbacks)
