# MERFISH Data
See [MERFISH poster PDF](../visualization/figures/MERFISH_HPC_Pipeline_Cowan_RCSymposium2024_poster.pdf) for a comprehensive overview and pilot data results of the computational workflow and **pilot study results** using a custom-designed 338-gene MERFISH panel.

- See poster code for [segmentation](01_vpt-segmentation) and [analysis](02_analysis/multi-sample/).

![MERFISH Poster Preview](../visualization/figures/MERFISH_HPC_Pipeline_Cowan_RCSymposium2024_poster.png)

## Code sections

- `01_vpt-segmentation`: This includes a series of scripts for re-segmenting our MERFISH data using a cellpose segmentation implementation of the Vizgen Post-Processing Tool (VPT). Both of the following subdirectories use the same general VPT workflow, with differences in the number of different segmentation applied to a single experiment.
    - `01_single-algorithm` Use these scripts as a workflow if you are interested in segmenting a single or small number of MERFISH datasets. Simple workflow; follow the scripts in order, `01_`, then `02_`, then, `03_`, then `04_`.
    - `02_hyperparameter-tuning` This section includes scripts I used when segmenting a single dataset many different ways to optimize segmentation results and understand how slight changes affect downstream clustering and computational analysis. Use these scripts as a template if you are testing the results of many different segmentation results. I use some bash scripting to automate this process. 
    
- `first-run-july-2022`: This section includes a Python script I used to explore BIG's first MERFISH dataset, generated in collaboration with Vizgen. This includes a very basic scripting workflow and the upstream data utilizes watershed segmentation. I highly recommend resegmenting MERFISH data with cellpose 2.0 because the default watershed segmentation supported by MERLIN has not performed as well in our hands for mouse brain.

- `02_analysis`: This includes the full computational workflow used to generate the pilot data reported in the linked poster.

## Designing a MERFISH Panel

We use the FPKM (abundance) values from bulk-sequencing designing MERFISH panels for highly-inflamed brain tissue. This is because MERFISH is a highly sensitive technique that is susceptible to optical crowding and compromised data collection. Becuase the brain is an immune-privileged site under homeostasis, the FPKM abundance of immune transcripts varies substantially between control and infected states, for instance. 

Here are guidelines from Vizgen for bulk-seq FPKM considerations
- 500 gene panel: FPKM must be under 8k
- 960-gene panel: FPKM limit is 15k

## Wet lab troubleshooting

- Make your best effort to maintain an RNase-free workstation
- To prevent tissue from shattering during the freezing process, we found that freezing in OCT made a very big dfiference
- During our initial runs we noticed significant inter-run variability. In conversations with Vizgen they have advised sectioning all samples to be included in an experiment once, store the samples in ethanol for up to 1 month, and proceeding with successive data collection runs.

## Considerations for data quality
In my discussions with Vizgen support, they have recommend >30k transcripts per FOV as indication of a successful MERFISH experiment. Lower transcript count can correspond to degraded RNA, and a good indication of this is a poor transcript density across FOV. Counterintuitively, I learned that **optical crowding** can be diagnosed by a reduction in transcripts per FOV relative to surrounding FOV.

![Figure 1](../visualization/figures/merfish-spatial-scatter.png)
