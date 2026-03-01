from pyspeedy import Speedy
from pyspeedy.callbacks import DiagnosticCheck, XarrayExporter
from pyspeedy.config import load_config

config = load_config()
run_config = config.run

model = Speedy(
    start_date=run_config.start_date,
    end_date=run_config.end_date,
    physics_mode="held_suarez",
)

# Held-Suarez uses generated flat lower-boundary data by default when no BC file is provided.
model.set_bc()

callbacks = [
    DiagnosticCheck(interval=run_config.diag_interval),
    XarrayExporter(
        interval=run_config.history_interval,
        output_dir=run_config.output_dir,
        variables=run_config.output_vars,
        verbose=run_config.verbose_output,
    ),
]

print(callbacks[1].variables)
model.run(callbacks=callbacks)
