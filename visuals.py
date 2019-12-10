import plotly.graph_objects as go
import numpy as np

# Create figure
fig = go.Figure()

data = open('data_18', 'r')
X = []
Y = []
c = 0
for line in data:
    if(c % 3 == 0):
        X.append([ line.split()[a] for a in range(0, len(line.split()), 3)])
        Y.append([ line.split()[a+1] for a in range(0, len(line.split()), 3)])
    c+=1
print (c)
for x in X:
    x.append(-30)
    x.append(30)

for y in Y:
    y.append(-30)
    y.append(30)

e_num = len(X)
# Add traces, one for each slider step
for epoch in np.arange(e_num):
    fig.add_trace(
        go.Scatter(
            visible=False,
            mode='markers',
            marker = dict(size = 3),
            name=str(epoch),
            x=X[epoch],
            y=Y[epoch]
        ))

# Make 10th trace visible
fig.data[10].visible = True

# Create and add slider
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="restyle",
        args=["visible", [False] * len(fig.data)],
    )
    step["args"][1][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=10,
    currentvalue={"prefix": "Epoch: "},
    pad={"t": 50},
    steps=steps
)]

fig.update_layout(
    sliders=sliders
)

fig.show()
