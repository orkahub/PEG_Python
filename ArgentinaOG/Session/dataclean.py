import pandas as pd

df_md = pd.read_csv('metadata.csv')
df_prod = pd.read_csv('production.csv')
df_hf = pd.read_csv('hf.csv')


dict_meta = {
    'idpozo': 'idwell',
    'sigla': 'wellname',
    'coordenadax': 'longitude',
    'coordenaday': 'latitude',
    'tipo_reservorio': 'reservoir_type',
    'cuenca': 'basin',
    'areapermisoconcesion': 'concession',
    'cota': 'elevation',
    'profundidad': 'depth',
    'fecha_data': 'dates'
}

dict_prod = {
    'idempresa': 'idcompany',
    'sigla': 'wellname',
    'anio': 'year',
    'mes': 'month',
    'idpozo': 'idwell',
    'prod_pet': 'oilprod',
    'prod_gas': 'gasprod',
    'prod_agua': 'waterprod',
    'iny_agua': 'waterinj',
    'iny_gas': 'gasinj',
    'iny_co2': 'co2inj',
    'tef': 'effective_time',
    'tipoextraccion': 'extraction_type',
    'empresa': 'company',
    'profundidad': 'depth',
    'formacion': 'formation',
    'cuenca': 'basin',
    'coordenadax': 'longitude',
    'coordenaday': 'latitude',
    'tipo_de_recurso': 'resource_type',
    'clasificacion': 'classification',
    'fecha_data': 'dates'
}

dict_hf = {
    'idpozo': 'idwell',
    'sigla': 'wellname',
    'tipo_reservorio': 'reservoir_type',
    'longitud_rama_horizontal_m': 'horizontal_length',
    'cantidad_fracturas': 'amount_fractures',
    'arena_bombeada_nacional_tn': 'national_sand',
    'arena_bombeada_importada_tn': 'imported_sand',
    'agua_inyectada_m3': 'injected_water',
    'co2_inyectado_m3': 'injected_co2',
    'presion_maxima_psi': 'max_pressure',
    'potencia_equipos_fractura_hp': 'equipment_power_hp',
    'fecha_data': 'dates'
}


def datacleaning(df_name, dict_name, filename):
    df_name.rename(columns=dict_name, inplace=True)
    kept_col_meta = [dict_name[x] for x in dict_name]
    df_name = df_name[kept_col_meta]
    df_name.to_csv(filename)

dict_list = [dict_meta, dict_prod, dict_hf]
df_list = [df_md, df_prod, df_hf]
filenames = ['metadata_eng.csv', 'production_eng.csv', 'hf_eng.csv']

for i, item in enumerate(dict_list):
    datacleaning(df_list[i], item, filenames[i])


print('')