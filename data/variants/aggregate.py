import pandas as pd
from datetime import date

start_date = date(2021, 1, 6)
end_date = date(2021, 3, 14)


def read_rows(_date: date) -> pd.DataFrame:
    df: pd.DataFrame = pd.read_csv(f'{_date.isoformat()}.csv').dropna()
    banned_states = ['Total']
    df = df[~df.State.isin(banned_states)]
    df['State'].replace('District of Columbia', 'DC')
    cols = df.columns.values.tolist()
    if 'Rate' in cols:
        return df.rename(columns={'State': 'state', 'Rate': 'cases'})[['state', 'cases']]
    if cols == ['State', 'Cases']:
        return df.rename(columns={'State': 'state', 'Cases': 'cases'})
    elif 'filter' in cols:
        return df[df['filter'] == 'Variant B.1.1.7'][['State', 'B.1.1.7 Variant ']].rename(columns={'State': 'state', 'B.1.1.7 Variant ': 'cases'})
    elif 'B.1.1.7 Variant ' in cols:
        return df[['State', 'B.1.1.7 Variant ']].rename(columns={'State': 'state', 'B.1.1.7 Variant ': 'cases'})


gdf = pd.DataFrame()
for _date in pd.date_range(start_date, end_date):
    try:
        df = read_rows(_date.date())
        if df is None:
            print(f'warning: {_date.date().isoformat()}')
            continue
        df['date'] = _date.date()
        gdf = gdf.append(df)
    except FileNotFoundError:
        print(f'omitted {_date.date().isoformat()}')

gdf.to_csv('variants.csv', index=False)
