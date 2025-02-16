#inaturalist observations
import pandas as pd
from pyinaturalist import get_observations, Observation

def get_inaturalist_observations(geometry) -> pd.DataFrame:
    response = get_observations(
    taxon_id=75400,       # species ID for Amaranthus tuberculatus
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
inat_obs_df = get_inaturalist_observations(us_bbox)
print(len(inat_obs_df))
