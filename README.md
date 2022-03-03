# Bayes project - Schools: ranking school examination results using multivariate hierarcical models

This project is conducted as part of the Bayesian Statistics and Hierarchical Models course by [Mathieu RIBATET](http://mribatet.perso.math.cnrs.fr/teaching.html#BAYES).

## Usage

### R

The R version of this project can be found in the directory ``R version``.

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
  poetry install
  poetry shell
  ```

Then, to run the MCMC algorithm you have to use the following command:
```sh
python main.py <infile>
```

where ``<infile>`` must be replaced by the path of the Json file e.g. ``data/data.json``.

## Authors

- Charles MIRANDA
- Armand ROUYRRE
- Aymard SAHNOUNE
- Vincent SEVESTRE


