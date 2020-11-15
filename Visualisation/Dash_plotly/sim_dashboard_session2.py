import dash
import dash_html_components as html
import dash_core_components as dcc
import pandas as pd
from dash.dependencies import Output, Input
import plotly.graph_objs as go
import numpy as np
import sqlalchemy as db
import json


qty_dict = {'PROD_gor': 'PROD - Gas Oil Ratio',
            'PROD_cumm_prod_gas': 'PROD - Cumulative Gas Production',
            'PROD_pw': 'PROD - Bottom Hole Pressure',
            'PROD_q_oilp': 'PROD - Oil Production Rate',
            'GINJ_cumm_inj_gas': 'GINJ - Cumulative Gas Injection',
            'GINJ_pw': 'GINJ - Bottom Hope Pressure',
            'GINJ_q_gasi': 'GINJ - Gas Injection Rate',
}
file = r'Sim Training sensitivities3.xlsx'
sheet_name = 'SPE1_PHYSICS'
file = 'Sim Training sensitivities3.xlsx'
df_input = pd.read_excel(file, sheet_name=sheet_name)

drop_down_options1 = [{'label': value, 'value': key} for key, value in qty_dict.items()]

input_list_filtered = ['CASE', 'struct_dipping', 'rock_compressibility', 'scal', 'Ginj', 'mob_ratio']
#input_list_filtered = list(df_input.columns)
drop_down_options2 = [{'label': value, 'value': value} for value in input_list_filtered]


# Get list of simulation cases
case_names = df_input['CASE'].tolist()

#Import Simulation results
engine = db.create_engine('postgresql://postgres:postgres@localhost:5432/SimResultsAdolfo')
exercise_name2 = "'" + 'Exercise 2' + "'"
metadata = db.MetaData()
df_master = pd.read_sql("select definition_table_id, name, description from master3 where description = 'Exercise 2'",
                       con=engine)
case_names = df_master['name'].tolist()
drop_down_options3 = [{'label': value, 'value': value} for value in case_names]

external_stylesheets = ['https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css']


app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    html.H1('Reservoir Simulation Dashboard'),
    html.Div(id='selected-data', style={'display': 'none'}),
    html.Div([
        html.Div([
            dcc.Dropdown(
                id='dd1',
                options=drop_down_options1,
                value='PROD_gor'
            ),
            dcc.Dropdown(
                id='dd2',
                options=drop_down_options3,
                value=case_names[0],
                multi=True,
            ),
        ], className='col-md-6'),
        html.Div([
            dcc.Dropdown(
                id='dd3',
                options=drop_down_options2,
                value='struct_dipping'
            ),
            dcc.Dropdown(
                id='dd4',
                options=drop_down_options2,
                value='rock_compressibility',
            ),

        ], className='col-md-6'),
    ],className='row'),
    html.Div([
        html.Div([
            dcc.Graph(
                id='graph1',
            )
        ],className='col-md-6'),
        html.Div([
            dcc.Graph(
                id='graph2',
            )

        ], className='col-md-6')
    ], className='row'),


])

def get_sim_names(selectedData):
    sim_names = []
    if selectedData == 'null':
        sim_names = None
    else:
        #selectedData = json.loads(selectedData)
        for item in selectedData['points']:
            sim_names.append("CASE_REV2_EX2_" + str(item['pointIndex']))
    return sim_names

@app.callback(Output('graph1', 'figure'),
             [Input('dd1', 'value'),
              Input('dd2', 'value')])

def update_plot(qty_plot, value2):
    df_temp = df_master[df_master['name'].isin(value2)]
    ressim_id = df_temp['definition_table_id'].tolist()
    all_ressim_id = "','".join(ressim_id)

    table_query = "select definition_table_id, quantity, values from third_table3 where definition_table_id in (" + "'" + all_ressim_id  + "'" + ") and quantity in ('time','" + qty_plot + "')"
    print(table_query)
    df_sim_results = pd.read_sql(
        table_query,
        con=engine)

    traces = []
    for i, case_id in enumerate(ressim_id):
        df_sim_results_temp = df_sim_results[df_sim_results['definition_table_id']==case_id]
        time_data = df_sim_results_temp[df_sim_results_temp['quantity']=='time']['values'].values
        qty_data = df_sim_results_temp[df_sim_results_temp['quantity']==qty_plot]['values'].values
        time_array = time_data[0].replace("[", "")
        time_array = time_array.replace("]", "")
        time_array = time_array.split(" ")
        # exec("qty_array = " + qty_array[0])
        qty_array = qty_data[0].replace("[ ", "")
        qty_array = qty_array.replace("[", "")
        qty_array = qty_array.replace(" ]", "")
        qty_array = qty_array.replace("]", "")
        qty_array = qty_array.split(" ")
        trace = go.Scatter(
            x=time_array,
            y=qty_array,
            name=value2[i],
            mode='lines'
        )
        traces.append(trace)

    figure = go.Figure(data=traces)
    return figure


@app.callback(Output('graph2', 'figure'),
             [Input('dd3', 'value'),
              Input('dd4', 'value')])
def update_second_plot(in1, in2):
    trace = go.Scatter(
        x=df_input[in1].values,
        y=df_input[in2].values,
        mode='markers'
    )
    traces = [trace]
    figure = go.Figure(data=traces)
    return figure


@app.callback(Output('selected-data', 'children'),
             [Input('graph2', 'selectedData')])
def store_selected_data(selectedData):
    sim_names = get_sim_names(selectedData)
    return sim_names


@app.callback(Output('dd2', 'value'),
             [Input('selected-data', 'children')])
def update_list_simcases(selectedCases):
    return selectedCases

if __name__ == '__main__':
    app.run_server(debug=True, port=8051)





