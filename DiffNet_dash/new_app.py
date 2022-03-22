import dash
import dash_cytoscape as cyto
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
import dash_table as dt
import plotly.express as px
from dash.dependencies import Input, Output, State
import pandas as pd
import sqlite3
import json
import random
import re
from dash.exceptions import PreventUpdate
import numpy as np
# from app import app

####################################
#### IMPORT database and query #####
####################################

data_path = './assets/'
meta = pd.read_csv(data_path + 'highly_filtered_features_edges_count.txt',sep='\t')
networks = pd.read_csv(data_path + 'highly-filtered-hits-12kEdges-133features-PCCge248-dPCCge100.txt',sep='\t').set_index('Condition')
networks.insert(0, 'Condition', list(networks.index))
layouts = json.loads(open(data_path + 'layouts.json','r').read())
cell_lines = pd.read_csv(data_path + 'features-bool-deDuped-2918feats-808cells.txt',sep='\t',index_col=0)

gene_profiles_path = './assets/gene_profiles.db'
source = 'ERBB2'
target = 'ERBB2'

with sqlite3.connect(gene_profiles_path) as gene_conn:
    profile_slice = pd.read_sql_query('SELECT * FROM profiles WHERE gene="{source}" OR gene="{target}"'.format(source=source,target=target),gene_conn).set_index('gene')
    profile_slice = profile_slice['bf'].apply(lambda x: list(json.loads(x).items())).explode()
    paired_profile = pd.DataFrame(profile_slice.tolist(),columns = ['depmap_id','bf'],index = profile_slice.index).reset_index().set_index(['gene','depmap_id']).unstack().dropna(axis=1).T
if gene_conn:
    gene_conn.close()


############  Layout  ##############
# 1. Title bar                     #
# 2. Left column:				   #	
#	 i)  	sidebar				   #
#	 ii) 	interactive network    #
#	 iii) 	table switch button	   #
#	 iv)	lab site link          #
# 3. Right column:                 #
#	 i)		 network table         #
#	 ii)     scatter plot          #
# 4. Link to Hart lab 			   #
####################################


external_stylesheets = ['./style.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = 'DiffNet'


####################################
####### Left column contents #######
####################################

gene_options = [{"label":g,"value":g} for g in np.unique(networks[['Gene1','Gene2']].values.flatten())]
gene_dropdown = dcc.Dropdown(
    options=gene_options,
    multi=True,
    id='gene_dropdown',
    style={'width':'30vw'}
    )
gene_selector = html.Div(
    children = [html.H3('Search for context-gene interactions',style={'margin':'0'}),gene_dropdown],
    style = {'width':'30vw','margin':'auto','paddingTop':'2vh'}
)
data_options = [
    {"label":'Tissue/tumor','value':'tissue/tumor'},
    {"label":'Genomic lesions','value':'genomic_lesions'}
    ]
data_dropdown = dcc.Dropdown(
    options = data_options,
    id = 'data_dropdown',
    value='genomic_lesions',
    style = {'width':'15vw','display':'inline-block'}
)

gene_bool = (meta['Condition'].str.contains('GOF') | meta['Condition'].str.contains('LOF'))
network_options = meta[gene_bool]['Condition'].sort_values().apply(lambda x: {"label":' '.join(x.split('_')),'value':x}).tolist()
network_dropdown = dcc.Dropdown(
    options = network_options,
    id = 'network_dropdown',
    style = {'width':'15vw','display':'inline-block'}
)
network_selector_header = html.Div(
    children = [
        html.H3('Functional context type',style={'margin':'0','width':'15vw','display':'inline-block'}),
        html.H3('Specific context',style={'margin':'0','width':'15vw','display':'inline-block'})
    ]
)
network_selector = html.Div(
    children = [network_selector_header,data_dropdown,network_dropdown],
    style = {'width':'30vw','margin':'auto','paddingTop':'2vh'}
)

min_range = dcc.Input(
    id="min_range",
    type="number",
    min=1,
    max=meta['Num_Edges'].max(),
    step=1,
    value=1,
    style={'width':'10vw','display':'inline-block','marginLeft':'-2.5vw'}
)
max_range = dcc.Input(
    id="max_range",
    type="number",
    min=1,
    max=meta['Num_Edges'].max(),
    step=1,
    value=meta['Num_Edges'].max(),
    style={'width':'10vw','display':'inline-block','marginLeft':'5vw'}
)
range_header = html.Div(
    children = [
        html.H3('Min. number of edges',style={'margin':'0','width':'10vw','display':'inline-block','marginLeft':'-2.5vw'}),
        html.H3('Max. number of edges',style={'margin':'0','width':'10vw','display':'inline-block','marginLeft':'5vw'})
    ]
)
network_range = html.Div(
    children = [range_header,min_range,max_range],
    style = {'width':'30vw','marginLeft':'10vw','paddingTop':'2vh'}
)


details_header = html.H3('Network details',style={'margin':'0','marginLeft':'5vw','display':'block','width':'20vw'})
network_details = html.Div(
    children = ['text'],
    style = {'height':'25vh','width':'34vw','backgroundColor':'white','margin':'auto'}
)
details = html.Div(
    children = [details_header,network_details],
    style = {'height':'25vh','width':'40vw','paddingTop':'2vh'}

)

sidemenu = html.Div(
    children = [gene_selector,network_selector],
    className = 'sidemenu',
    style = {'height':'25vh','width':'40vw','marginTop':'2vh','background-color': '#aaaaaa'} #'#d3d3d3'
)

left_panel = html.Div(
    children = [sidemenu],
    style = {'height':'25vh','width':'40vw','marginLeft':'5vw','display':'inline-block', 'top':0}

)

bottom_width = 40 #same as left panel width

refresh_network_button = html.Button(['Refresh network'],id='refresh',style={'position':'absolute','zIndex':1},n_clicks=0)
graph = cyto.Cytoscape(
        id='cytoscape',
        layout={'name': 'preset'},
        style={
            'width': '{left_width}vw'.format(left_width=bottom_width),
            'height': '30vw',
            'outlineStyle':'solid','outlineWidth':'medium','outlineColor':'#aaaaaa',
            },
        elements=[],
        zoom=.1,
        minZoom=.1,
        # panningEnabled = False,
        maxZoom=5,
        responsive=False,
        stylesheet = [
        {
            'selector': '[weight < 0]',
            'style': {
                'line-style':'dashed'
            }
        },
        # Group selectors
        {
            'selector': 'node',
            'style': {
                'content': 'data(label)',
                'font-size':20
            }
        },
        ]
    )

interactive_title = html.Div([
					  
					        html.H2('Interactive network',
					        style={
					            'width':'{left_width}vw'.format(left_width=bottom_width),
					            'height':'5vh',
					            'textAlign':'center',
					            'paddingTop':'2vh'
					            })					        
        ])


switch = html.Div(children = [daq.BooleanSwitch(
  id='table_switch',
  label='Click on for a tabular view of network',
  labelPosition='top',
  on=False,
  style = {'display':'block', 'textAlign':'center'}
)])

lab_link = html.Div( id='link_to_Hartlab', children=[
    				dcc.Link( 'Hart Lab', href='https://www.hart-lab.org')
    				], style={'width':'40vw','textAlign':'center',
    						  'font-size': '18px'}
    			)

bottom_panel = html.Div(
    children = [interactive_title, switch, refresh_network_button, graph, lab_link],
    className='bottom_panel',
    id='bottom_panel',
    style={
        'width':'40vw','height':'35vh',
        'display':'block',
        'position':'absolute',
        #'left':'50vw',
        'marginLeft':'5vw',
        'top':'35vh'
        }
)

########################################################################################################################################


####################################
####### Right column contents ######
####################################

########## scatter plot ############

right_width = 40
right_height = 40
scatter = px.scatter(paired_profile,x = source,y = target)
t = scatter.update_layout(
    plot_bgcolor='white',
    yaxis={'zeroline':True,'zerolinecolor':'black','zerolinewidth':2,'showgrid':False,'linecolor':'black','linewidth':2},
    xaxis={'zeroline':True,'zerolinecolor':'black','zerolinewidth':2,'showgrid':False,'linecolor':'black','linewidth':2}
)
figure = dcc.Graph(id='figure',figure=scatter,style={'display':'none'})
plot = html.Div(children=figure,id = 'plot_container',style={'width':'{right_width}vw'.format(right_width=right_width), 
															 'display':'flex', 'justify-content': 'center', 'align-items': 'center'})
##############################



############# table #############

table_width = 40

net_table = html.Div(
    children = [
        dt.DataTable(
            id='meta_table',
            columns=[],
            data=None,
            filter_action='native',
            sort_action='native',
            page_size=50,
            style_table={'width':'{table_width}vw'.format(table_width=table_width - 2),'display':'none'}
        )],style={
        	'top':0,
            'height':None, # to avoid the scrollbar
            'overflowY':'auto',
            'clear':'both',
            'display':'block'}
)


right_panel = html.Div(
    children = [plot,net_table],
    className='right_panel',
    id='right_panel',
    style={
    	'top':0,
        'width':'{right_width}vw'.format(right_width = right_width + 1),'height':'{right_height}vh'.format(right_height=right_height),
        'position':'absolute',
        'left':'{left_space}vw'.format(left_space = 10+bottom_width + 2),
        'display':'inline-block',
        'paddingTop':'12vh'
        #'overflowX':'hidden'
        }
)


################################################################################################################################################

app.layout = html.Div([

####################################
############## Title ###############
####################################

			html.H1('Dynamic functional interactions', style={'textAlign':'center'}),


####################################
####### Left column contents #######
####################################


			html.Div([
				    left_panel,
				    bottom_panel
				    ]),

####################################
####### Right column contents ######
####################################


			html.Div(
				    children = [right_panel],
				    style={
				        'width':'40vw',
				        'height':'40vh',
				        'overflowX':'hidden'#,
				        },
				    className='panels'
				    ),

###################################
######## Additional items #########
###################################

			html.Div(
				    [
				        html.P(id='node_click',children=[None],style={'display':'none'}),
				        html.P(id='edge_source_click',children=[None],style={'display':'none'}),
				        html.P(id='edge_target_click',children=[None],style={'display':'none'})
				    ],
				    style={'display':'block'}
				    )


		], style={'display':'inline-block',
				  'width':'100vw',
				  'height':'100vh'})


########################################################################################################################################

#####################################
############ Callbacks ##############
#####################################


@app.callback(
    Output("gene_dropdown", "options"),
    [Input("gene_dropdown", "search_value")],
)
def update_options(search_value):
    if not search_value:
        opt = gene_options.copy()
    else:
        opt = [o for o in gene_options if re.match(search_value, o["label"], re.IGNORECASE)]
        opt.extend([o for o in gene_options if o not in opt and search_value in o["label"]])
    return opt

@app.callback(
    Output("network_dropdown", "options"),
    [Input("network_dropdown", "search_value"),Input('gene_dropdown','value'),Input('data_dropdown','value')],
)
def update_options(search_value,genes,dataset):

    if dataset == 'tissue/tumor':
        default_options = meta[~gene_bool]['Condition'].sort_values().apply(lambda x: {"label":' '.join(x.split('_')),'value':x}).tolist()
    else:
        default_options = meta[gene_bool]['Condition'].sort_values().apply(lambda x: {"label":' '.join(x.split('_')),'value':x}).tolist()

    if not search_value:
        options = default_options.copy()
    else:
        options = [o for o in default_options if re.match(search_value, o["label"], re.IGNORECASE)]
        options.extend([o for o in default_options if o not in options and search_value in o["label"]])

    if genes != None and len(genes) > 0:
        networks_gene_search = networks['Gene1'].isin(genes) | networks['Gene2'].isin(genes)
        conditions = networks.loc[networks_gene_search,'Condition'].drop_duplicates().tolist()
        options = [o for o in options if o['value'] in conditions]

    if not search_value and genes == None:
        options = default_options.copy()

    return options



@app.callback(
    [Output('cytoscape','elements'),Output('meta_table','columns'),Output('meta_table','data')], #Output('interactive_link','href')],
    [Input('network_dropdown','value'),Input('refresh', 'n_clicks')],
    prevent_initial_call=True
    )
def update_elements(key, button):
    if key in layouts.keys():
        elements = json.loads(layouts[key])
        for e in elements:
            e['locked'] = True
    else:
        elements = []

    if key == '':
        raise PreventUpdate
    else:
        network = networks.loc[[key]].copy()
        network.round(3)
        network['Gene pairs'] = network['Gene1'] + '-' +  network['Gene2']
        header = ['Gene pairs'] + list(set(network.columns) - {'Gene pairs'})

    return([elements, [{"name": i, "id": i} for i in network[header].columns],network.round(4).to_dict('records')])#,'/apps/dash_net/' + key])

def set_scaling(N):
    scale = pd.DataFrame({5:200,25:500,50:750,200:1000,500:2500,4000:5000}.items(),columns = ['size','scale'])
    min_scale = scale[scale['size'].ge(N)]['size'].min()
    index = scale[scale['size'].eq(min_scale)].index[0]
    return(scale.loc[index,'scale'])

@app.callback(
    Output('node_click', 'children'),
    [Input('cytoscape', 'tapNodeData')]
    )
def displayTapNodeData(data):

    if data:
        return(data['label'])
    else:
        return(None)

@app.callback(
    [Output('edge_source_click', 'children'),Output('edge_target_click', 'children')],
    Input('cytoscape', 'tapEdgeData')
    )
def displayTapEdgeData(data):
    if data:
        return(data['source'],data['target'])
    else:
        return(None,None)

@app.callback(
    Output('plot_container','children'),
    [Input('edge_source_click','children'),Input('edge_target_click','children')],
    [Input('network_dropdown','value')]
)
def update_figure(source,target,key):
    network_key = key
    if source == None or target == None:
        return(None)
    else:
        with sqlite3.connect(gene_profiles_path) as conn:
            c = conn.cursor()
            profile_slice = pd.read_sql_query('SELECT * FROM profiles WHERE gene="{source}" OR gene="{target}"'.format(source=source,target=target),conn).set_index('gene')
            profile_slice = profile_slice['bf'].apply(lambda x: list(json.loads(x).items())).explode()
            paired_profile = pd.DataFrame(profile_slice.tolist(),columns = ['depmap_id','bf'],index = profile_slice.index).reset_index().set_index(['gene','depmap_id']).unstack().dropna(axis=1).T
        if conn:
            conn.close()

        sampled_data = set(cell_lines.index[cell_lines[network_key]].tolist()).intersection(set(paired_profile.index.levels[1]))

        paired_profile['group'] = 'out'
        paired_profile.loc['bf'].loc[sampled_data,'group'] = 'in'
        if 'LOF' in network_key:
            color_map = {'in':'#1f77b4','out':'grey'}
        elif 'GOF' in network_key:
            color_map = {'in':'#ff7f0e','out':'grey'}
        else:
            color_map = {'in':'#2ca02c','out':'grey'}
        scatter = px.scatter(paired_profile,x = source,y = target,color='group',trendline='ols',color_discrete_map=color_map)
        scatter.update_layout(
            #plot_bgcolor='white',
            paper_bgcolor='rgba(0,0,0,0)',
    		plot_bgcolor='rgba(0,0,0,0)',
            yaxis={'zeroline':True,'zerolinecolor':'black','zerolinewidth':2,'showgrid':False,'linecolor':'black','linewidth':2},
            xaxis={'zeroline':True,'zerolinecolor':'black','zerolinewidth':2,'showgrid':False,'linecolor':'black','linewidth':2}
        )
        plot = [dcc.Graph(id='figure',figure=scatter,style={'display':'inline-block','width':'30vw','height':'30vw'})] # set the scatter plot to a square
        return(plot)

@app.callback(
    Output('meta_table','style_table'),
    [Input('table_switch','on')]
)
def toggle_table(switch):
    if switch:
        style={
            'height':None, # to avoid the table scrollbar
            'overflowY':'hidden',
            'clear':'both',
            'display':'inline-block'}
    else:
        style={
            'height':None,
            'overflowY':'hidden',
            'clear':'both',
            'display':'none'}
    return(style)

if __name__ == '__main__':
	app.run_server(debug=True)





	
