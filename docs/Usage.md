# Usage

### Docker Installation
#### Whole Genome Variant Calling

For whole genome variant calling, or calling variants on several chromosomes, use `NanoCaller_WGS.py` to call variants. Assuming all your input files are in a folder `YOUR_INPUT_DIR`, and you want to use `YOUR_OUTPUT_DIR` to store the results. 
```
VERSION=0.3.0
docker run -itd \
-v 'YOUR_INPUT_DIR':'/input/' \
-v 'YOUR_WORKING_DIR':'/output/' \
umahsn/nanocaller:${VERSION} \
python NanoCaller_WGS.py \
-bam /input/YOUR_BAM \
-ref /input/YOUR_REF \
-prefix PREFIX \
-o /output \
-cpu NUMBER_OF_CPUS_TO_USE 
```

#### Single Chromosome Variant Calling

For calling variants on single chromosomes, use `NanoCaller.py` to call variants. Assuming all your input files are in a folder `YOUR_INPUT_DIR`, and you want to use `YOUR_OUTPUT_DIR` to store the results. 
```
VERSION=0.3.0
docker run -itd \
-v 'YOUR_INPUT_DIR':'/input/' \
-v 'YOUR_WORKING_DIR':'/output/' \
umahsn/nanocaller:${VERSION} \
python NanoCaller.py \
-chrom CHROMOSOME \
-bam /input/YOUR_BAM \
-ref /input/YOUR_REF \
-prefix PREFIX \
-o /output \
-cpu NUMBER_OF_CPUS_TO_USE 
```
