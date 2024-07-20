library(Seurat)
print('Start adding new annotation')
print(levels(Idents(seuratObjFinal)))

seuratObjFinal     <- RenameIdents(seuratObjFinal,
                                   ## ---
                                   `ST1` = 'Stromal',
                                   `ST2` = 'Stromal',
                                   `ST3` = 'Stromal',
                                   `ST4` = 'Stromal',
                                   `ST5` = 'Stromal',
                                   `SM` = 'Smooth mucle',
                                   `CE` = 'Celiatic epthelial',
                                   `SE` = 'Secretory epthelial',
                                   `P/V1` = 'Pericytes',
                                   `P/V2` = 'Pericytes',
                                   `P/V3` = 'Pericytes',
                                   `EN1` = 'Endothelial',
                                   `EN2` = 'Endothelial',
                                   `EN3` = 'Endothelial',
                                   `EN4` = 'Endothelial',
                                   `LE` = 'Lymphatic',
                                   `T/NK1` = 'T/NK',
                                   `T/NK2` = 'T/NK',
                                   `T/NK3` = 'T/NK',
                                   `MP` = 'Macrophage',
                                   `MA` = 'Mast'
                                   ## ---
)
