# Physio Analysis

- `code/physio.py` contains the functions necessary to currently
  process the GSR data in Bike_ATrain and BikeExtend
- `code/explore_physio.ipynb`, is a notebook that contains the "playground"
  for thinking about how to organize the analysis.
- `environment.yml` specifies the requisite python packages necessary to run the code.
- `data` contains example test data for evaluating the functions.
- `output` contains example output from the script.

## Example call

From the root directory of this repository:

```bash
./code/physio.py data output/testOut
```

To get help type: `./code/physio.py --help`

Example Docker call:

```bash
docker run -it --rm -v ${PWD}:/home/neuro/physioAnalysis -p 8445:8080 jdkent/physio-analys:dev
```

