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
	CONDA_DIR=$(shell conda info --base)
	ENV_DIR=$(CONDA_DIR)/envs/$(ENVIRONMENT_NAME)
    CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
endif

environment_create:
ifeq (True,$(HAS_CONDA))
ifneq ("$(wildcard $(ENV_DIR))","")
	@echo ">>> Found $(ENVIRONMENT_NAME) environment in $(ENV_DIR). Skipping installation..."
else
	@echo ">>> Detected conda, but '$(ENVIRONMENT_NAME)' environment is missing. Installing ..."
	conda env create --file $(ENVIRONMENT_FILE)
	@echo "> Activate environment: '$(ENVIRONMENT_NAME)'"
	$(CONDA_ACTIVATE) $(ENVIRONMENT_NAME)
	@echo "> Show environment list"
	conda env list
endif
else
	@echo "> Install conda first"
endif

environment_remove:
	@echo "> Uninstalling the '$(ENVIRONMENT_NAME)' environment"
	conda remove --yes --name $(ENVIRONMENT_NAME) --all
	@echo "> Show environment list"
	conda env list