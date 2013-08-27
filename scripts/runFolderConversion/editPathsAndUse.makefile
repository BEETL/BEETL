BEETL_CONVERT=/home/ljanin/install/BEETL-release/bin/beetl-convert
INPUT_DIR=/illumina/scratch/optimusprime/EagleSimulations/20130702_Sergii_30x30x/EAGLE/EAGLE_normal/RunFolder/Data/Intensities/BaseCalls

LANES=$(patsubst $(INPUT_DIR)/L00%, %, $(wildcard $(INPUT_DIR)/L00*))
CYCLES=$(patsubst $(INPUT_DIR)/L001/C%.1, %, $(wildcard $(INPUT_DIR)/L001/C*.1))
TILES=$(patsubst $(INPUT_DIR)/L001/C1.1/s_1_%.bcl, %, $(wildcard $(INPUT_DIR)/L001/C1.1/s_1_*.bcl))

$(warning LANES:$(LANES))
$(warning TILES:$(TILES))
$(warning CYCLES:$(CYCLES))


all: $(foreach cycle,$(CYCLES),cyc.$(cycle))

cyc.%:
	$(foreach lane,$(LANES), $(foreach tile,$(TILES), \
	  $(BEETL_CONVERT) --input-format=bcl --output-format=cyc -i $(INPUT_DIR)/L00$(lane)/C$*.1/s_$(lane)_$(tile).bcl -o cyc_$(lane)_$(tile).$* ; \
	  cat cyc_$(lane)_$(tile).$* >> $@ ; \
	  rm -f cyc_$(lane)_$(tile).$* ; \
	  rm -f cyc_$(lane)_$(tile).$*.qual ; \
	)) \
	echo "ok $@"


# following doesn't seem to work yet
sge: $(foreach cycle,$(CYCLES),sge.cyc.$(cycle).done)

sge.cyc.%.done:
	$(MAKE) -n cyc.$* | sed -e 's/make/#make/g' > sge.cyc.$*.script \
	&& sync \
	&& sync \
	&& qsub -sync y -cwd -v PATH sge.cyc.$*.script \
	&& touch $@
