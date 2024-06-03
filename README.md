# Bayesian-based fire models

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Obtaining Driving Data](#obtaining-driving-data)
- [Configuration Settings](#configuration-settings)
- [Model Setup](#model-setup)
- [Running the Models](#running-the-models)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)



## Introduction
This repository contains a series of Bayesian-based fire models. These models follow a common inference and sampling workflow to predict fire behaviour based on various driving data. . At the moment, we have two models that target burned area:
*	ConFire, has been tested and is working as expected.
*	FLAME, is also operational but is located in a separate repository. It needs to undergo testing in our current environment. You can find it at https://github.com/malu-barbosa/FLAME.
  
Additionally, we are in the process of working on INFERNO, which is derived from JULES-ES-INFERNO ([Mangeon et al. 2016]( https://gmd.copernicus.org/articles/9/2685/2016/)) and working on a  fireMIP ensemble optimisation scheme. We are also extending these two models to cover other aspects of fire regimes beyond burned area.

## Features
- Bayesian inference for fire modelling
- Flexible configuration settings
- Different fire models for relating drivers to burnt area
- Scalable to different data sets and scenarios
- Detailed logging and result output

## Installation
Follow these steps to install the necessary software and dependencies:

1. Clone the repository:
    ```bash
    git clone [https://github.com/yourusername/your-repo.git](https://github.com/douglask3/Bayesian_fire_models.git)
    ```
2. Navigate into the project directory:
    ```bash
    cd Bayesian_fire_models
    ```
3. Install the required packages:
    
    ```bash
    conda env create -f environment.yml
    conda activate bayesian-fire-models
    ```

## Obtaining Driving Data
The driving data required for the models are listed [here](https://github.com/douglask3/Bayesian_fire_models/blob/main/README/Datasets). Once you have downloaded the dataset you want, them save it to `data/data/driving_data/` directory within the project folder, or update the path in the corrisponding configution file (see below) to its location.

## Configuration Settings
Configuration settings can be found in the `config.yaml` file. Here are the key settings you may need to adjust:

- `parameter1`: Description
- `parameter2`: Description
- `parameter3`: Description

## Model Setup
Each model follows the same basic setup process. Here are the steps:

1. Navigate to the model directory:
    ```bash
    cd models/model_name
    ```
2. Adjust the configuration file if necessary.
3. Prepare the environment by running:
    ```bash
    setup_script.sh
    ```

## Running the Models
To run a model, execute the following command:

```bash
python run_model.py --config config.yaml



#  Fogo local Analisado pela Máxima Entropia (FLAME)

# Controls on Fire (ConFire)

# Note for State of Wildfire report


# FLAME -  Fogo local Analisado pela Máxima Entropia

## Introduction
FLAME (Fogo local Analisado pela Máxima Entropia) is a Bayesian inference implementation of a maximum entropy fire model specifically tailored to simulating fires in heterogeneous territories like Brazil. 

## Model overview

### Maximum Entropy concept

### Relationship curves

### Model Optimization

## Evaluation overview

### Posterior analysis

### Jackknife

### Response curves

### NME

## Projections

### Attribution

### Out of sample 

## Data processing

### Basic datasets provided

### Notes on constraining area


#### *Political areas*
When opening data, setting ''Country'' or ''Continent'' will constrain the extent to that country or continent, and mask areas outside of it. Uses Natural Earth. If you define a country, it one look at the continent. Use None is you don't want any. Continent options are:
* 'South America'
* 'Oceania'
* 'Europe'
* 'Afria'
* 'North America'
* 'Asia'

#### *Ecoregions*
''ecoregions'' is a numeric list (i.e [3, 7, 8]) where numbers pick Olson bomes and mask out everywhere else. If you  pick more than one, it returns a map of all of them.
* **None** Return all areas.
* **1** Tropical and subtropical moist broadleaf forests
* **2** Tropical and subtropical dry broadleaf forests
* **3** Tropical and suptropical coniferous forests
* **4** Temperate broadleaf and mixed forests
* **5** Temperate Coniferous Forest
* **6** Boreal forests / Taiga
* **7** Tropical and subtropical grasslands, savannas and shrublands
* **8** Temperate grasslands, savannas and shrublands
* **9** Flooded grasslands and savannas
* **10** Montane grasslands and shrublands
* **11** Tundra
* **12** Mediterranean Forests, woodlands and scrubs
* **13** Deserts and xeric shrublands
* **14** Mangroves

#### *Brazillian legal Biomes*

''Biomes'' is a numeric list where numbers pick Brazilian biomes and mask out everywhere else. If you pick more than one, it returns a map of all of them.

* **1** Amazonia
* **2** Caatinga
* **3** Cerrado
* **4** Atlantic Forest
* **5** Pampa
* **6** Pantanal

#### *GFED Regions*

#### *AR6 regions*

#### *To year range*

#### *To months of year*

