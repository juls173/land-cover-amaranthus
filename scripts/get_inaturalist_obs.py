#inaturalist observations
import pandas as pd
from pyinaturalist import get_observations, Observation

def get_inaturalist_observations(geometry, taxon_id: int) -> pd.DataFrame:
    """
    Retrieves iNaturalist observations of Amaranthus tuberculatus within a specified geographic region.

    This function queries the iNaturalist API for observations of Amaranthus tuberculatus (taxon ID: 75400) 
    that have open location data. It filters results based on the provided geographic geometry 
    and returns the data as a pandas DataFrame, and then saves the results to a csv file.

    Args:
        geometry: A geospatial geometry object defining the area of interest.

    Returns:
        pd.DataFrame: A DataFrame containing observation dates and locations.
    """
    response = get_observations(
    taxon_id=taxon_id,       # species ID for Amaranthus tuberculatus
    geoprivacy='open',    # only include observations with open location data
    geoframe=geometry,  
    page='all'
)
    observations = Observation.from_json_list(response)

    obs_data = []
    for obs in observations:
        obs_data.append({
        'date': obs.observed_on,
        'location': obs.location,
    })

    obs_df = pd.DataFrame(obs_data)
    obs_df.to_csv('data/inat_amaranthus_tuberculatus_observations.csv', index=False)
    return obs_df

us_bbox = (24.6, -124.8, 49.0, -66.9)
inat_obs_df = get_inaturalist_observations(us_bbox, taxon_id=75400)
print(len(inat_obs_df))
