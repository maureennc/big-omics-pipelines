# 7/30/24 MC

## Examples of correct and incorrect subsetting + how to verify

##################################################################################################################

# Tumor - A vs. B
## Incorrect
Tumor_A_B <- subset(target_data, select = phenoData(target_data)[["class2"]] == c("Tumor_A", "Tumor_B"))

## Correct
Tumor_A_B <- target_data[, pData(target_data)[["class2"]] %in% c("Tumor_A", "Tumor_B")]

## Verification
table(pData(target_data)[["class2"]], pData(target_data)[["slide"]])
table(pData(Tumor_A_B)[["class2"]], pData(Tumor_A_B)[["slide"]])

##################################################################################################################

# Normal tissue - A vs. B
## Incorrect
Normal_A_B <- subset(target_data, select = phenoData(target_data)[["class2"]] == c("Normal_A", "Normal_B"))

## Correct
Normal_A_B <- target_data[, pData(target_data)[["class2"]] %in% c("Normal_A", "Normal_B")]

## Verification
table(pData(target_data)[["class2"]], pData(target_data)[["slide"]])
table(pData(Normal_A_B)[["class2"]], pData(Normal_A_B)[["slide"]])

##################################################################################################################

# Normal_B vs. Tumor_B
## Incorrect
Normal_B_Tumor_B <- subset(target_data, select = phenoData(target_data)[["class2"]] == c("Tumor_B", "Normal_B"))

## Correct
Normal_B_Tumor_B <- target_data[, pData(target_data)[["class2"]] %in% c("Normal_B", "Tumor_B")]

## Verification
table(pData(target_data)[["class2"]], pData(target_data)[["slide"]])
table(pData(Normal_B_Tumor_B)[["class2"]], pData(Normal_B_Tumor_B)[["slide"]])

##################################################################################################################
