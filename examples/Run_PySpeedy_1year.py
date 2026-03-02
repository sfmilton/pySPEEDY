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


config = load_config()
run_config = config.run
model_config = config.model

start_date = run_config.start_date
end_date = _add_calendar_year(start_date)
summary_interval = 30 * model_config.nsteps

# Create an instance of the speedy model using the defaults in pyspeedy/data/model_config.yml,
# but extend the run window to one calendar year. `end_date` is exclusive.
model = Speedy(
    start_date=start_date,
    end_date=end_date,
)

# Initialize the boundary conditions from the default bundled dataset.
model.set_bc()

callbacks = [
    DiagnosticCheck(interval=run_config.diag_interval),
    RuntimeSummary(
        interval=summary_interval,
        verbose=run_config.verbose_output,
    ),
    XarrayExporter(
        interval=run_config.history_interval,
        output_dir=run_config.output_dir,
        variables=run_config.output_vars,
        verbose=run_config.verbose_output,
    ),
    DailyDiagnosticsExporter(
        output_dir=run_config.output_dir,
        verbose=run_config.verbose_output,
    ),
]
history_exporter = callbacks[2]

print(f"Running from {start_date:%Y-%m-%d} to {end_date:%Y-%m-%d} (end date exclusive).")
print(f"Runtime summaries every {summary_interval // model_config.nsteps} days.")
print(history_exporter.variables)

model.run(callbacks=callbacks)
