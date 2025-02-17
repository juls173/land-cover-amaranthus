from pygbif import occurrences
import pandas as pd

def get_amaranthus_occurences(bounding_box: str) -> pd.DataFrame:
    """
    Retrieves occurrences of Amaranthus tuberculatus within a specified bounding box 
    from the GBIF database.

    Args:
        bounding_box (str): A WKT (Well-Known Text) polygon defining the search area.

    Returns:
        pd.DataFrame: A DataFrame containing occurrence records with longitude, latitude, 
                      event date, and basis of record.
    """
    data = occurrences.search(
    taxon_key = 8577467,
    geometry=bounding_box,
    hasCoordinate=True,
    limit=1000,
    basisOfRecord="PRESERVED_SPECIMEN",
    year="2000,2023", # range of years, can narrow down
)

    df = pd.DataFrame(data['results'])
    df = df[['decimalLongitude', 'decimalLatitude', 'eventDate', 'basisOfRecord']]

    df.to_csv('data\\amaranthus_occurrences.csv', index=False)
    return df

il_bbox = "POLYGON((-91.513 36.970, -87.495 36.970, -87.495 42.508, -91.513 42.508, -91.513 36.970))"
all_bbox = "POLYGON(())"
df = get_amaranthus_occurences(il_bbox)
print(df.iloc[:5])