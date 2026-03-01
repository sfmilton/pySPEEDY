from pyspeedy import Speedy
from pyspeedy.callbacks import DiagnosticCheck, XarrayExporter
from pyspeedy.config import load_config

config = load_config()
run_config = config.run

# Create an instance of the speedy model using the defaults in pyspeedy/data/model_config.yml.
model = Speedy(
    start_date=run_config.start_date,
    end_date=run_config.end_date,
)
# At this point, the model state is "empty".

# To initialize the model, we need to define its boundary conditions first.
# This function will set the default boundary conditions derived from the ERA reanalysis.
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

# Print the names of output variables that will be saved.
# Note that the variables shown next are in the grid space (not the spectral space)
print(callbacks[1].variables)

# Run the model
model.run(callbacks=callbacks)
# After the model is run, the model state will keep the last values of the last integration step.
