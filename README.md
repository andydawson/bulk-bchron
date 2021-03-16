# bulk-bchron
Workflow to fit age-depth models for North American fossil pollen records.

To run the workflow that fits age-depth models for the set of North American fossil pollen cores in Neotoma Paleoecology Database that still have dates in Radiocarbon years BP using Bchron, run `run_bchron.r`. This code will try to fit an age-depth for each site that has a folder (with the datasetid as the folder name) in the Cores folder.

A PDF of a figure showing the estimated age-depth model as well as the chronological constraints is saved in the datasetid folder (`datasetid_bchron.pdf`). A compilation of the PDFs for all of the sites for which age-depth models were fit is saved in the main working directory (`bchron_plots_vX.X.pdf`). A report indicating the success or failure of the age-depth model fitting process is compiled in `bchron_report_vX.X.csv`.

From the age-depth model, posterior samples of the ages are drawn for each pollen depth (`datasetid_bchron_samples.csv`) as well as for the depths of the chronological constraints (`datasetid_bchron_geo_samples.csv`).

