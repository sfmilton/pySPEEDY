"""
This scripts is used to regenerate the fixtures.
It should only be used when some parts of the model change and the reference simulations
needs to be recomputed.
"""

from datetime import datetime

from pyspeedy import Speedy
from pyspeedy.callbacks import DiagnosticCheck, XarrayExporter
from pyspeedy.config import load_config

config = load_config()
run_config = config.run

model = Speedy(
    start_date=datetime(1982, 1, 1),
    end_date=datetime(1982, 1, 4),
)
model.set_bc()

callbacks = [
    DiagnosticCheck(interval=run_config.diag_interval),
    XarrayExporter(
        interval=run_config.history_interval,
        output_dir=run_config.output_dir,
        verbose=run_config.verbose_output,
    ),
]

model.run(callbacks=callbacks)
