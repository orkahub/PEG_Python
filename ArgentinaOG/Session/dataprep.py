import pandas as pd

url_metada = 'http://datos.minem.gob.ar/dataset/c846e79c-026c-4040-897f-1ad3543b407c/resource/cbfa4d79-ffb3-4096-bab5-eb0dde9a8385/download/listado-de-pozos-cargados-por-empresas-operadoras.csv'
url_wellprod = 'http://datos.minem.gob.ar/dataset/c846e79c-026c-4040-897f-1ad3543b407c/resource/b5b58cdc-9e07-41f9-b392-fb9ec68b0725/download/produccin-de-pozos-de-gas-y-petrleo-no-convencional.csv'
url_hf = 'http://datos.minem.gob.ar/dataset/71fa2e84-0316-4a1b-af68-7f35e41f58d7/resource/2280ad92-6ed3-403e-a095-50139863ab0d/download/datos-de-fractura-de-pozos-de-hidrocarburos-adjunto-iv-actualizacin-diaria.csv'
url_list = [url_metada, url_wellprod, url_hf]
filenames = ['metadata.csv', 'production.csv', 'hf.csv']


def dataIO(url, filename):
    df = pd.read_csv(url)
    df.to_csv(filename)


for i, item in enumerate(url_list):
    dataIO(item, filenames[i])




