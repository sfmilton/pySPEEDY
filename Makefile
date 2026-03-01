# Main makefile for Speedy.f90
PYTHON ?= python

.PHONY:clean
all: 
	$(PYTHON) scripts/render_model_config.py
	$(MAKE) -C speedy.f90
	@cp speedy.f90/libspeedy.a ./

.PHONY:clean
test:	
	pytest -s -v pyspeedy

.PHONY:clean
clean:   
	$(MAKE) -C speedy.f90 clean
	@rm -f pyspeedy/*.so
	@rm libspeedy.a
