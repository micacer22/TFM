import pandas as pd

# Crear la tabla de condiciones
samples = ["AI1", "AI2", "AI3", "AM1", "AM2", "AM3", 
           "CI1", "CI2", "CI3", "CM1", "CM2", "CM3",
           "SI1", "SI2", "SI3", "SM1", "SM2", "SM3",
           "WI1", "WI2", "WI3", "WM1", "WM2", "WM3"]

conditions = ["atx1_infected"] * 3 + ["atx1_mock"] * 3 + \
             ["clf_infected"] * 3 + ["clf_mock"] * 3 + \
             ["sdg8_infected"] * 3 + ["sdg8_mock"] * 3 + \
             ["WT_infected"] * 3 + ["WT_mock"] * 3

# Crear DataFrame
metadata_df = pd.DataFrame({"sample": samples, "condition": conditions})

# Guardar como archivo CSV
metadata_path = "/Users/miriamcaballerocervero/Desktop/metadata_conditions.csv"
metadata_df.to_csv(metadata_path, index=False)

# Devolver la ruta del archivo guardado
metadata_path

