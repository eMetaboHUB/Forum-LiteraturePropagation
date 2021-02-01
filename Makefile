.PHONY: environment

# Oneshell means I can run multiple lines in a recipe in the same shell, so I don't have to
# chain commands together with semicolon
.ONESHELL:
# Need to specify bash in order for conda activate to work.
SHELL=/bin/bash
ENVIRONMENT_FILE = env/propagation.yml
ENVIRONMENT_NAME = propagation

ifeq (,$(shell which conda))
    HAS_CONDA=False
else
    HAS_CONDA=True
    CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
endif

environment_create:
ifeq (True,$(HAS_CONDA))
	@echo "> Detected conda, creating conda environment '$(ENVIRONMENT_NAME)' (Overwrite if already exists)"
	conda env create --force --file $(ENVIRONMENT_FILE)
	@echo "> Activate environment: '$(ENVIRONMENT_NAME)'"
	$(CONDA_ACTIVATE) $(ENVIRONMENT_NAME)
	@echo "> Show environment list"
	conda env list
else
	@echo "> Install conda first"
endif

environment_remove:
	@echo "> Uninstalling the '$(ENVIRONMENT_NAME)' environment"
	conda remove --yes --name $(ENVIRONMENT_NAME) --all
	@echo "> Show environment list"
	conda env list