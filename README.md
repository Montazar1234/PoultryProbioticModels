# Metabolic Modeling of Fungal Probiotics in the Poultry Gut Microbiome


![image](https://github.com/Montazar1234/PoultryProbioticModels/assets/22119956/a2025388-4eca-4367-a4fc-a8c7953e518b)


## Repository Structure

- `Data/`: Contains metagenomic datasets, as well as manifest files used for simulation.
- `micom_script.py`: Python script for running MiCOM simulations.
- `data/`: Directory containing input data files.
- Metabolic agora models (https://github.com/Gibbons-Lab/isb_course_2023)
- CarveFungi metabolic models (https://zenodo.org/records/7413265)
- `Results/`: Directory containing output files.

## Dependencies

The code in this repository relies on the following dependencies:

- MiCOM [(version 0.33.2)](https://micom-dev.github.io/micom/)
- CarveFungi ((https://github.com/SandraCastilloPriego/CarveFungi)
- AutoPACMEN https://github.com/klamt-lab/autopacmen

Please refer to the `requirements.txt` file for the complete list of dependencies and their versions.

## Usage

1. Clone this repository: `git clone https://github.com/your-username/your-repository.git`
2. Install the required dependencies: `pip install -r requirements.txt`
3. Prepare the input data files and place them in the `data/` directory.
4. Run the Jupyter Notebook `metabolic_modeling_poultry.ipynb` to reproduce the analysis and generate the results.
5. Alternatively, run the MiCOM simulation script: `python micom_script.py`

## Data

The input data files used in this study are located in the `data/` directory. The files include:

- `aggregated_relative_abundances.ods`: Metagenomic data of the poultry gut microbiome.
- `manifest.csv`: Manifest file containing information about the fungal strains.
- `diet_composition.csv`: Composition of the poultry diets used in the simulations.


## Contact

For questions or inquiries, please contact

- Montazar Al-Nijir (Man69@bath.ac.uk)

## Acknowledgments

We would like to thank the University of Bath's Research Computing Group for their support in this work.
We would also like to thank Sandra Castillo for her help and advice during the use of CarveFungi and Paula Jouhten for her valuable advice and support. 
Finally this repository contains code adapted from the notebook "MiCOM 2023" by Gibbons Lab, available at: [https://colab.research.google.com/github/Gibbons-Lab/isb_course_2023/blob/main/micom_2023.ipynb]
