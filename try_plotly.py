# (*) To communicate with Plotly's server, sign in with credentials file
import plotly.plotly as py

py.sign_in('vasiliy.triandafilidi', 'h5uejpbuf4')

# (*) Useful Python/Plotly tools
import plotly.tools as tls

# (*) Graph objects to piece together plots
from plotly.graph_objs import *

import numpy as np  # (*) numpy for math functions and arrays
import matplotlib.pyplot as plt


# Package all mpl plotting commands inside one function
def plot_mpl_fig():

    # Make two time arrays
    t1 = np.arange(0.0, 2.0, 0.1)
    t2 = np.arange(0.0, 2.0, 0.01)

    # N.B. .plot() returns a list of lines.
    # The "l1, = plot" usage extracts the first element of the list
    # into l1 using tuple unpacking.
    # So, l1 is a Line2D instance, not a sequence of lines
    l1, = plt.plot(t2, np.exp(-t2), label='decaying exp.')
    l2, = plt.plot(t2, np.sin(2 * np.pi * t2), '--go', label='sin')
    l3, = plt.plot(t1, np.log(1 + t1), '.', label='log')
    l4, = plt.plot(t2, np.exp(-t2) * np.sin(2 * np.pi * t2), 'rs-.',label='exp')

    # Add axis labels and title
    plt.xlabel('time')
    plt.ylabel('volts')
    plt.title('Damped oscillation')
    return (l1, l2, l3, l4)  # return line objects (for legend, later)

# Plot it!



# (l1, l2, l3, l4) = plot_mpl_fig()
# plt.legend((l2, l4), ('oscillatory', 'damped'), 'upper right')

# mpl_fig1 = plt.gcf()

# py_fig1_ss = tls.mpl_to_plotly(mpl_fig1, strip_style=True)

# py.plot(py_fig1_ss, filename='s6_damped_oscillation-default-style2',
#         auto_open=False)

# py_fig1.strip_style()

# py.plot(py_fig1, filename='s6_damped_oscillation-default-style3',
#         auto_open=False)

# # N.B. get matplotlib figure object and assign a variable to it

# Convert mpl fig obj to plotly fig obj, resize to plotly's default
# plot_mpl_fig()
# (l1, l2, l3, l4) = plot_mpl_fig()
# plt.legend((l2, l4), ('oscillatory', 'damped'), 'upper right')
# get current figure, retrieve the information on the figure that has been just created
#alternative mpl_fig1 = plt.figure()
# mpl_fig2 = plt.gcf()
# # mpl_fig2 = plt.figure()
# py_fig2 = tls.mpl_to_plotly(mpl_fig2, resize=True)

# # Delete misplaced legend annotations
# py_fig2['layout'].pop('annotations', None)

# # Add legend, place it at the top right corner of the plot
# py_fig2['layout'].update(
#     showlegend=True,
#     legend=Legend(
#         x=1.05,
#         y=1
#     )
# )

# # Send updated figure object to Plotly, show result in notebook
# # py.iplot(py_fig2, filename='s6_damped_oscillation-legend2')

# py.plot_mpl(py_fig2, strip_style=True,
#              filename='s6_damped_oscillation-default-style')



#use the function described above to plot everything
plot_mpl_fig()
# get current figure, retrieve the information on the figure that has been just created
#alternative mpl_fig1 = plt.figure()
mpl_fig1 = plt.gcf()

#now we created a plotly object and we are printing the information about it
py_fig1 = tls.mpl_to_plotly(mpl_fig1, verbose=True)

# print information
print(py_fig1.to_string())
print py_fig1
for i in ['autosize', 'width', 'height']:
    print i, py_fig1['layout'][i]

py_fig1['layout'].pop('annotations', None)

# Add legend, place it at the top right corner of the plot
py_fig1['layout'].update(
    showlegend=True,
    legend=Legend(
        x=1.05,
        y=1
    )
)

py_fig1['data'][0].update(name='oscillatory0')
py_fig1['data'][1].update(name='oscillatory1')
py_fig1['data'][2].update(name='oscillatory2')
py_fig1['data'][3].update(name='damped3')

# Do not include the remaining traces in legend
py_fig1['data'][0].update(showlegend=True)
py_fig1['data'][1].update(showlegend=False)
py_fig1['data'][2].update(showlegend=False)
py_fig1['data'][2].update(showlegend=False)

py.plot_mpl(mpl_fig1, strip_style=True,
             filename='s6_damped_oscillation-default-style')


# send the matplotlib figure to online
# plot_url = py.plot_mpl(mpl_fig1)

# py.iplot_mpl(mpl_fig1, filename='s6_damped_oscillation')
# py.iplot_mpl(mpl_fig1, strip_style=True,
#              filename='s6_damped_oscillation-default-style')
