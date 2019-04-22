# LivinOnTheEdge
Lake connectivity for aquatic and semi-aquatic species;
not: https://en.wikipedia.org/wiki/Livin'_on_the_Edge

# ConceptualModel
GIS and powerpoint files used to make conceptual model

# Data
Input and output datasets used in analysis; also contains metadata

# ItFigures
Output figures for manuscript

# Rcode
Scripts for data analysis and output figures

functions: contains custom function for lake and wetland statistics within dispersal buffer around focal lake
Hydro_Terrestrial_LakeConn_Indices.R: bulk of analysis for manuscript; calculates aquatic and semi-aquatic connectivity indices for each lake using principal component analysis (using data from other scripts). Generates most of figures in manuscript.
id_lakeconn_lakes.R: used to identify lakes hydrologically connected to Great Lakes
LakesWetlands_inDispersalBuffers.R: executes custom function for lake and wetland statistics within dispersal buffer around focal lake
Michigan_LAGOS_conn_metrics.R: extract freshwater connectivity metrics for Michigan lakes from LAGOS
Michigan_up_down_lakeconn_metrics.R: calculate additional freshwater connectivity metrics not in LAGOS
MichiganLakePatchStats.R: calculate "patch" statistics for Michigan lakes, treating lakes like terrestrial habitat patches
PercentProtected_IWS_LakeBuffers.R: calculate percent protected for lake watersheds and dispseral buffers (using tabulate area output from ArcGIS)
