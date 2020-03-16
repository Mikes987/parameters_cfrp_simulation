# Results of ROFit
This folder summarizes the results and thus files created by ```ROFit.py```. After the script is finished, one can observe:
- Diagrams
- Chart of operating temperatures and Ramberg-Osgood Parameters
- Logfile that documents how the data was treated.

## Diagrams
Each png file shows two diagrams.
- Left: All single real shear tests as well as their mean curce.
- Right: Mean curve and its associated Ramberg-Osgood fit.

## Logfile
In order to create a sufficient mean curve, two requirements were necessary:
- All datasets for a specific temperature shall have the same number of xy datapoints
- These datapoints shall be as narrow as possible

Both requirements were not guaranteed by the equipment that was used, so an adjustment had to be realized. That was done by using the following steps:
- For each operating temperature, the dataset with the smallest number of datapoints is set to be the reference number.
- First, all other datasets with their xy-datapoints beginning at index ```len(reference)``` are examined. All datapoints with x-values bigger than the highest x-value of the reference will be deleted.
- If the number of datapoints is still uneven loop through each line and remove the minimum in each line as long as the number of datapoints becomces even. It is simply assumed that the minimum in each line is far smaller than the other values.

Each loop is documented in the logfile.
