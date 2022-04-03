### Code to add controls to list b/c if in main code, rows will be added multiple times ###

add_euro = data.frame ('Pre-EuroAmerican settlement horizon', '0', 'NA', '0')
write.table(add_euro, file = "chroncontrol_types-edited.csv", sep = ",",
            append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE)

add_euro2 = data.frame ('European settlement horizon', '0'', 'NA', '0')
write.table(add_euro, file = "chroncontrol_types-edited.csv", sep = ",",
            append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE)

picea_decline = data.frame ('Picea decline', '0', 'NA', '0')
write.table(picea_decline, file = "chroncontrol_types-edited.csv", sep = ",",
            append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE)

picea_decline2 = data.frame ('IDW-2d Picea decline', '0', 'NA', '0')
write.table(picea_decline2, file = "chroncontrol_types-edited.csv", sep = ",",
            append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE)


