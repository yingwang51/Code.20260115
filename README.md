## README：Data and R code for “Cultivation interacts with life form to mediate the effect of native range size on plant naturalization”

#### 

### data

**1. 20250330.native_plant.matched.cultivated_plant_DO.csv：** 

**Description**: This dataset constitutes the core data for analysing the global naturalisation patterns of native Chinese plants. It integrates species naturalization information, cultivation information, life forms, native range size, and cultivation  both within and outside China.

#Taxonomic Identifiers

- `accepted_plant_name_id`: The unique id for the accepted name

- `taxon_name`: The accepted name (Latin name) of the species

- `taxon_authors`: The authors who named the species

#Life forms & Cultivation status & naturalization status

- `life.form.integrated`: Integrated life form of the species 

- `cultivation.lin2019`: Cultivation status determined according to* Catalogue of Cultivated Plants in China*（yes or no）

- `nat.extent`: Naturalization extent (number of regions where plants are naturalized)

- `nat.incidence`: Naturalization incidence (likelihood of becoming naturalized)。

#native range size & cultivation extent

- `native.global.tdwg3`: Native range size, (i.e., number of regions where native)

- `planting.China.tdwg3`: Cultivation within China, measured as the number of TDWG level-3 regions within China where the species is recorded as cultivated.

- `WorldCuP.n.tdwg3`: Cultivation outside China, measured as the number of TDWG level-3 regions outside China where the species is recorded as cultivated.

**2. NativeChinesePlantsCultivationExtentWorldCuP.csv**

**Description**: This dataset quantifies the cultivation extent of Chinese native plants outside of China. Out of the 25,808 native taxa analyzed, 9,670 are identified as being cultivated internationally.

- `taxon_name`: The accepted name (Latin name) of the species

- `taxon_authors`: The authors who named the species

- `accepted_plant_name_id`: The unique id for the accepted name

- `number_records`: The number of unique WorldCuP lists in which each taxon is present

- `number_tdwg3`: The total number of unique TDWG level-3 regions outside of China where the taxon is recorded as cultivated according to WorldCuP

**3. References_WorldCuPSources.csv***

**Description**: This file contains the comprehensive list of 365 bibliographic references and data sources used to compile the cultivation extent information for Chinese native plants within the World Checklist of Cultivated Plants (WorldCuP) database.

### code

**1.Rcode for Figure_2.R**

#Visualisation of BSEM analysis

**2.Rcode for Figure_S1.R**

#The distribution of naturalized plant species that are native to China across TDWG level-3 regions, separately for each life form (annual herbaceous, perennial herbaceous, and woody)

**3.Rcode for Figure_S3.R**

#Linearized relationships from Hurdle Models

**4.Rcode for Figure_S4 and Table_S2.R**

#BSEM analysis and posterior distributions of pathway-specific effects estimated from BSEMs

**5.Rcode for Figure_S6.R**

#correlation analysis (Native Range Size, Cultivation within/outside China)

**6.Rcode for Figures_1,S2 and Table_S1.R**

#Hurdle model analyse






