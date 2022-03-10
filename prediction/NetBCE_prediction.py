# %%
import numpy as np
import os,sys
import argparse
import numpy as np
from collections import defaultdict
from keras.models import Model, load_model, model_from_json
from datetime import datetime
import plotly.graph_objects as go
import base64
from dominate.util import raw
from dominate.tags import *
from dominate import document
np.random.seed(12345)
os.environ['CUDA_VISIBLE_DEVICES']="-1"

# %%
#binding model
main_dir = "../models"
def import_model(main_dir):
    models = []
    json_f = open(main_dir + "/model.json", 'r')
    loaded_model_json = json_f.read()
    json_f.close()
    for i in range(5):
        loaded_model = model_from_json(loaded_model_json)
        loaded_model.load_weights((main_dir + "/model_"+str(i)+".h5"))
        models.append(loaded_model)  
    return models

def scoring(models, data):
    import numpy as np
    '''
    Use an ensemble of models to yield prediction scores for the data
    Args:
        1. models: Model ensemble
        2. data: Data to be predicted
    Return values:
        1. probas_: Prediction scores
    '''
    probas_ = [np.transpose(model.predict(data))[0] for model in models]
    probas_ = [np.mean(scores) for scores in zip(*probas_)]
    return probas_     

# %%
def read_fasta(fasta_file):
    try:
        fp = open(fasta_file)
    except IOError:
        exit()
    else:
        fp = open(fasta_file)
        lines = fp.readlines()
        fasta_dict = {} 
        gene_id = ""
        for line in lines:
            if line[0] == '>':
                if gene_id != "":
                    fasta_dict[gene_id] = seq
                seq = ""
                gene_id = line.strip() #  line.split('|')[1] all in > need to be id
            else:
                seq += line.strip()        
        fasta_dict[gene_id] = seq #last seq need to be record       
    return fasta_dict 

# %%
def sample_fasta_peptides(sequences, peptide_lengths):
    sample_peptides = {}
    data_dict = defaultdict(defaultdict)
    for (i, (name, sequence)) in enumerate(sequences.items()):
        if not isinstance(sequence, str):
            raise ValueError("Expected string, not %s (%s)" % (
                sequence, type(sequence)))
        for peptide_start in range(len(sequence) - min(peptide_lengths) + 1):
            for peptide_length in peptide_lengths:
                peptide = sequence[peptide_start: peptide_start + peptide_length]
                start_stop = '_'.join([str(peptide_start + 1), str(peptide_start + peptide_length)])
                if len(peptide) != peptide_length:
                    continue
                if name not in sample_peptides.keys() :
                    sample_peptides[name] = [peptide]
                else:
                    sample_peptides[name].append(peptide)
                if name not in data_dict.keys() or peptide not in data_dict[name]:
                    data_dict[name][peptide] = [start_stop]
                else:
                    data_dict[name][peptide].append(start_stop)

    return sample_peptides, data_dict 

# %%
def convertSampleToProbMatr(sampleSeq3DArr): #changed add one column for '1'
    """
    Convertd the raw data to probability matrix
    PARAMETER
    ---------
    sampleSeq3DArr: 3D numpy array
       X denoted the unknow amino acid.
    probMatr: Probability Matrix for Samples. Shape (nb_samples, 1, nb_length_of_sequence, nb_AA)
    """

    letterDict = {}
    letterDict["A"] = 0
    letterDict["C"] = 1
    letterDict["D"] = 2
    letterDict["E"] = 3
    letterDict["F"] = 4
    letterDict["G"] = 5
    letterDict["H"] = 6
    letterDict["I"] = 7
    letterDict["K"] = 8
    letterDict["L"] = 9
    letterDict["M"] = 10
    letterDict["N"] = 11
    letterDict["P"] = 12
    letterDict["Q"] = 13
    letterDict["R"] = 14
    letterDict["S"] = 15
    letterDict["T"] = 16
    letterDict["V"] = 17
    letterDict["W"] = 18
    letterDict["Y"] = 19
    letterDict["-"] =20
    AACategoryLen = 21 

    probMatr = np.zeros((len(sampleSeq3DArr), len(sampleSeq3DArr[0]), AACategoryLen))
    sampleNo = 0
    for sequence in sampleSeq3DArr:
        AANo = 0
        for AA in sequence:
            if not AA in letterDict:
                probMatr[sampleNo][AANo] = np.full((1,AACategoryLen), 1.0/AACategoryLen)
            else:
                index = letterDict[AA]
                probMatr[sampleNo][AANo][index] = 1
            AANo += 1
        sampleNo += 1
    del sampleSeq3DArr
    return probMatr

def convertSampleToProbMatr_2(sample_peptides):
    seq_matr_Dict = {}
    for seq in sample_peptides.keys():
        allele_data = sample_peptides[seq]
        all_seqs = []
        for ind, seq1 in enumerate(allele_data):
            my_len = len(seq1)
            if my_len < 24:
                add = 24 - my_len
                seq2 = seq1 + '-' * add
            all_seqs.append(seq2)
        seq_matr = convertSampleToProbMatr(all_seqs)
        if seq not in seq_matr_Dict.keys():
            seq_matr_Dict[seq] = seq_matr
    return seq_matr_Dict

# %%
def preiction(fastafile, output_folder, peptide_lengths):
    sequences = read_fasta(fastafile)
    sample_peptides, sample_peptides_position = sample_fasta_peptides(sequences, peptide_lengths)
    sample_peptides_matr = convertSampleToProbMatr_2(sample_peptides)
    model_ligand = import_model(main_dir)

    all_epi_scores = []
    for seq in sample_peptides_matr.keys():
        allele_data = sample_peptides_matr[seq]
        allele_data = np.array(allele_data)
        epi_scores = scoring(model_ligand, allele_data)
        all_peptide = sample_peptides[seq]
        for i in range(len(all_peptide)):
            pos = sample_peptides_position[seq][all_peptide[i]]
            all_epi_scores.append([seq,
                                all_peptide[i],
                                ';'.join(pos).replace('_', ':'),
                                str(epi_scores[i])])

    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    with open(output_folder + '/NetBCE_predictions.tsv', 'w') as f:
        for line in all_epi_scores:
            f.write('\t'.join(line) + '\n')    

    return all_epi_scores                             

# %%
def wrap_plotly_fig(fig: go.Figure, width: str = '100%', height: str = '100%'):
    if 'px' in width:
        fig = fig.to_html(include_plotlyjs=False, full_html=False, default_height=height, default_width=width)
        return div(raw(fig), style=f'width: {width}')
    else:
        fig = fig.to_html(include_plotlyjs=False, full_html=False, default_height=height, default_width='100%')
        return div(raw(fig), style=f'width: {width}')

def ploty_fig_to_image(fig: go.Figure, width: int = 360, height: int = 360):
    fig_data = fig.to_image(format='svg', width=width, height=height).decode()
    return img(src=f'data:image/svg+xml;base64,{fig_data}',
               className='img-fluid',
               style=f'width: 100%; height: auto')

def get_plotlyjs():
    fig = go.Figure()
    fig = fig.to_html(include_plotlyjs=True, full_html=False)
    plotlyjs = fig[fig.index("<script"):fig.rindex("<div id=")] + "</div></script>"
    return raw(plotlyjs)

def lab_logo():
    lab_logo = base64.b64encode(
        open(str('../logo/logo.png'), 'rb').read()).decode()
    return img(src=f'data:image/jpg;base64,{lab_logo}', className='img-fluid',
                style="max-width:100%; max-height:100%; margin-left: 10px;"
                        "margin-bottom: 8px")  # can add opacity: 50% to style if desired

def prediction_table(all_epi_scores, className=None):
    t = table(className=f'display nowrap',
                style="text-align: center",
                id='table_id_example')
    t.add(
        thead(
            tr(
                [
                    th(f'Protein identifier', style="padding: 5px"),
                    th(f'Candidate epitope', style="padding: 5px"),
                    th(f'Position (start:end)', style="padding: 5px"),
                    th('Predicion score', style="padding: 5px"),
                ]
            )
        )
    )
    tablebody = tbody()
    for sample in all_epi_scores:
        tablerow = tr()
        tablerow.add(td(sample[0], style='word-break: break-word'))
        tablerow.add(td(sample[1]))
        tablerow.add(td(sample[2]))
        tablerow.add(td(sample[3]))
        tablebody.add(tablerow)

    t.add(tablebody)
    return div(t, className=f'table-responsive {className}' if className else 'table-responsive')

def gen_prediction_histogram(all_epi_scores, className=None):

    n_peps_fig = go.Figure()
    all_length = [len(sample[1]) for sample in all_epi_scores]
    binders, counts = np.unique(all_length, return_counts=True)

    n_peps_fig.add_trace(go.Bar(x=binders, y=counts, marker=dict(color = binders,
                     colorscale='Spectral')))
    n_peps_fig.update_layout(margin=dict(l=20, r=20, t=20, b=20),
                                hovermode='x',
                                legend=dict(orientation="h",
                                            yanchor="bottom",
                                            y=1.02,
                                            xanchor="right",
                                            x=1,
                                            bgcolor="rgba(255, 255, 255, 0.8)"),
                                font_color='#212529'
                                )
    n_peps_fig.layout.title.xanchor = 'center'
    n_peps_fig.update_yaxes(title_text='Number of peptides')
    n_peps_fig.update_xaxes(title_text='Peptides length')
    n_peps_fig.update_xaxes(titlefont={'size': 16}, tickfont={'size': 14})
    n_peps_fig.update_yaxes(titlefont={'size': 16}, tickfont={'size': 14})
    # n_peps_fig.write_image(str(fig_dir / 'binding_histogram.pdf'), engine="kaleido")
    card = div(div(b('Epitope length distribution'), className='card-header'), className='card')
    card.add(div(raw(n_peps_fig.to_html(full_html=False, include_plotlyjs=False)), className='card-body'))
    return div(card, className=className)

def gen_prediction_boxplot(all_epi_scores, className=None):

    n_peps_fig = go.Figure()

    length_score_dict = {}
    for i, prediction in enumerate(all_epi_scores):
        pep_len = len(prediction[1])
        pep_score = float(prediction[-1])
        if pep_len not in length_score_dict.keys():
            length_score_dict[pep_len] = [pep_score]
        else:
            length_score_dict[pep_len].append(pep_score)

    binders = list(length_score_dict.keys())
    counts = list(length_score_dict.values())
    N = len(binders)
    colors = ['hsl('+str(h)+',50%'+',50%)' for h in np.linspace(0, 360, N)]

    for xd, yd, cls in zip(binders, counts, colors):
            n_peps_fig.add_trace(go.Box(
                y=yd,
                name=xd,
                boxpoints='all',
                jitter=0.5,
                whiskerwidth=0.2,
                fillcolor=cls,
                marker_size=2,
                line_width=1)
            )

    n_peps_fig.update_layout(
        margin=dict(l=20, r=20, t=20, b=20),
                    hovermode='x',
                    legend=dict(orientation="h",
                                yanchor="bottom",
                                y=1.02,
                                xanchor="right",
                                x=1,
                                bgcolor="rgba(255, 255, 255, 0.8)"),
                                font_color='#212529')

    n_peps_fig.layout.title.xanchor = 'center'
    n_peps_fig.update_yaxes(title_text='Score')
    n_peps_fig.update_xaxes(title_text='Peptides length')
    n_peps_fig.update_xaxes(titlefont={'size': 16}, tickfont={'size': 14})
    n_peps_fig.update_yaxes(titlefont={'size': 16}, tickfont={'size': 14})
    # n_peps_fig.write_image(str(fig_dir / 'binding_histogram.pdf'), engine="kaleido")
    card = div(div(b('Score distribution'), className='card-header'), className='card')
    card.add(div(raw(n_peps_fig.to_html(full_html=False, include_plotlyjs=False)), className='card-body'))
    return div(card, className=className)

# %%
def generate_html_report(all_epi_scores, peptide_lengths, output_folder):
    doc = document(title='NetBCE Report')
    with doc.head:
        link(rel="stylesheet", href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/css/bootstrap.min.css",
                integrity="sha384-9aIt2nRpC12Uk9gS9baDl411NQApFmC26EwAOH8WgZl5MYYxFfc+NcPb1dKGj7Sk",
                crossorigin="anonymous")
        link(rel="stylesheet", href="https://cdn.datatables.net/1.11.4/css/jquery.dataTables.min.css")
        link(rel="stylesheet", href="https://cdn.datatables.net/buttons/2.2.2/css/buttons.dataTables.min.css")
        script(src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js")
        script(src="https://maxcdn.bootstrapcdn.com/bootstrap/4.5.0/js/bootstrap.min.js")
        script(src="https://cdn.datatables.net/1.11.4/js/jquery.dataTables.min.js")
        script(src="https://cdn.datatables.net/buttons/2.2.2/js/dataTables.buttons.min.js")
        script(src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js")
        script(src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/pdfmake.min.js")
        script(src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/vfs_fonts.js")
        script(src="https://cdn.datatables.net/buttons/2.2.2/js/buttons.html5.min.js")
        script(src="https://cdn.datatables.net/buttons/2.2.2/js/buttons.print.min.js")

        body(onload="plots = document.getElementsByClassName('plotly-graph-div');"
                    "l = plots.length;"
                    "for (i=0; i < l; i++) {Plotly.relayout(plots[i], {autosize: true});}")
        script("$(document).ready(function(){$('.nav-tabs a').click(function(){$(this).tab('show');});"
                "$('.nav-tabs a').on('shown.bs.tab',"
                "function(){"
                "plots = document.getElementsByClassName('plotly-graph-div');"
                    "l = plots.length;"
                    "for (i=0; i < l; i++) {Plotly.update(plots[i]);}"
                "});});")
        # script("$(document).ready(function(){$('#table_id_example').DataTable();});")
        script("$(document).ready(function(){$('#table_id_example').DataTable({dom: 'Bfrtip',buttons: ['copy', 'csv', 'excel', 'pdf', 'print'], aaSorting: [[3, 'desc']]});});")

    with doc:
        get_plotlyjs()
        with div(id='layout', className='container', style='max-width: 1600px;'
                                                            'min-width: 1000px;'
                                                            'margin-top: 20px;'
                                                            'margin-bottom: 20px'):
            with div(className='row'):
                with div(className='col-12', style='display: flex; height: 60px'):
                    div([h1('N'), h3('et'), h1('BCE'),
                            h5(f' (v1)', style="white-space: pre"),
                            h1(' - Analysis report')],
                        style="background-color: #0c0c0c; padding: 5px; color: white;"
                                "border-radius: 6px; width: 100%; display: flex"),
                    lab_logo()
            hr()
            with div(className='row'):
                with div(className='col-10', style='margin: 0'):
                    h3('NetBCE enables accurate prediction of linear B-cell epitopes with interpretable deep neural network')
                    p('The identification of B-cell epitopes is of great value for the development of specific serodiagnostic assays and the optimization of medical therapy.Here, we present NetBCE, a python tool which uses a deep neural network to detect linear B-cell epitope regions on individual protein sequences. NetBCE exceeds all other currently used linear B-cell epitope prediction tools. Our software is shown to reliably predict linear B-cell epitopes of a given protein sequence, thus contributing to a significant reduction of laboratory experiments and costs required for the conventional approach.', style="font-size: 20px; padding: 5px;")
                    p([b('Developers: '), f'Haodong Xu, Zhongming Zhao @ CPH, UTHealth-Houston SBMI.'])
                    p([b('Lab website: '), a(f'https://www.uth.edu/bioinfo/.', href='https://www.uth.edu/bioinfo/')])
                    p([b('Date: '), f'{str(datetime.now().date())}'])
            hr()
            h3('Predicion results:')
            hr()
            with div(className='row', style="max-height:600px; overflow-y: scroll;"):
                with div(className='col-10', style='margin: 10px'):
                    prediction_table(all_epi_scores, className='col')

            hr()
            h3("Predicion statistics:")
            hr()
            with div(className='row'):
                if len(peptide_lengths) <= 12:
                    gen_prediction_histogram(all_epi_scores, className='col-6')
                    gen_prediction_boxplot(all_epi_scores, className='col-6')
                else:
                    gen_prediction_histogram(all_epi_scores, className='col-10')
                    gen_prediction_boxplot(all_epi_scores, className='col-10')
            hr()

    loc = f'{str("%s/report.html" % output_folder)}'
    with open(loc, 'w') as f:
        f.write(doc.render().replace("&lt;", "<"))        

# %%
def main():
    parser = argparse.ArgumentParser(description="progrom usage")
    parser.add_argument("-f", "--fasta", type=str, help="The input of antigen proteins")
    parser.add_argument("-l", "--length", type=int, default=[16,15,12,20], nargs='*', help="The length range of peptides")
    parser.add_argument("-o", "--out", type=str, help="prediction output path")
    args = parser.parse_args()
    fastafile = args.fasta
    peptide_lengths = args.length
    output_folder = args.out
    print('Job start')
    all_epi_scores = preiction(fastafile, output_folder, peptide_lengths)
    generate_html_report(all_epi_scores, peptide_lengths, output_folder)
    print('Job finish')
if __name__ == '__main__':
        main()

# %%
# python /collab2/hxu6/B_cell_epitope/python/NetBCE_prediction.py -f '/collab2/hxu6/B_cell_epitope/data/seq/test.fasta' -l 16 15 14 20 21 -o '/collab2/hxu6/B_cell_epitope/prediction'


