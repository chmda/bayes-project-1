# Bayes project - Schools: ranking school examination results using multivariate hierarcical models

This project is conducted as part of the Bayesian Statistics and Hierarchical Models course by [Mathieu RIBATET](http://mribatet.perso.math.cnrs.fr/teaching.html#BAYES).

:information_source: **If you want to read a more understandable code, please read the R version**.

## Usage

### R

The R version of this project can be found in the directory ``R version``.

Dependencies : 
  - ``matlib`` package

### Python

If the file ``data/data.json`` does not exist you need to run ``data/to_json.sh`` to convert the R data to a JSON file so that Python can access it easily.

- Python (no venv)

  First, you need to install the required packages using the following command:
  ```sh
  pip install -r requirements.txt
  ```
- Poetry
  
  If you use poetry, it's easy:
  ```sh
  poetry install --no-dev
  poetry shell
  ```

Then, to run the MCMC algorithm you have to use the following command:
```sh
python main.py <datafile> <initfile>
```

where ``<datafile>`` must be replaced by the path of the Data Json file e.g. ``data/data.json`` and ``<initfile>`` by the path of the Init Json File e.g. ``data/init1.json``.

## Results

After running the python version with the following command

```bash
python main.py data/data.json data/init1.json -n 10000 -b 1000
```

We obtained the following results

|  coefficient  |    mean   |    std   |   q2.5pc  |   median  |  q97.5pc |
|:--------:|:---------:|:--------:|:---------:|:---------:|:--------:|
|  beta[1] |  1.08e-03 | 6.87e-04 | -5.86e-04 |  1.08e-3  |  2.03e-3 |
|  beta[2] |  3.69e-01 | 1.34e-01 |  1.22e-01 |  3.69e-01 | 6.43e-01 |
|  beta[3] |  2.73e-01 | 8.28e-02 |  1.03e-01 |  2.73e-01 | 4.29e-01 |
|  beta[4] |  7.09e+00 | 3.05e+00 |  2.28e+00 |  7.09e+00 | 1.51e+01 |
|  beta[5] |  3.64e-01 | 3.81e-01 | -2.75e-01 |  3.64e-01 | 1.36e+00 |
|  beta[6] | -1.96e-01 | 6.68e-01 | -2.21e+00 | -1.96e-01 | 8.42e-01 |
|  beta[7] | -1.09e+00 | 1.56e+00 | -5.15e+00 | -1.09e+00 | 4.92e-01 |
|  beta[8] | -5.51e-01 | 3.35e-01 | -1.09e+00 | -5.51e-01 | 1.66e-01 |
| gamma[1] |  6.74e-03 | 1.01e+00 | -1.99e+00 |  6.74e-03 | 1.99e+00 |
| gamma[2] |  1.93e-03 | 9.92e-01 | -1.94e+00 |  1.93e-03 | 1.96e+00 |
| gamma[3] | -3.05e-02 | 9.91e-01 | -1.93e+00 |  -3.05e-2 | 1.88e+00 |
|    phi   | -8.19e-03 | 7.02e-03 | -1.57e-02 | -8.19e-03 | 9.37e-03 |
|   theta  | -3.85e-02 | 2.49e-01 | -8.28e-01 | -3.85e-02 | 2.21e-01 |
## Authors

- Charles MIRANDA
- Armand ROUYRRE
- Aymard SAHNOUNE
- Vincent SEVESTRE


