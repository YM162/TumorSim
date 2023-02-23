using DrWatson
@quickactivate "TumorSim"

using TumorSim
using Blink

w = Window()
body!(w,String(read("scripts/probando.html")),async=false)
loadjs!(w,"https://cdn.plot.ly/plotly-2.18.0.min.js")


js(w, Blink.JSString("""
                    var trace2 = {
                        x:[1,2,3], y: [1,2,3], z: [1,2,3],
                        mode: 'markers',
                        marker: {
                            color: 'rgb(127, 127, 127)',
                            size: 12,
                            symbol: 'circle',
                            line: {
                            color: 'rgb(204, 204, 204)',
                            width: 1},
                            opacity: 0.8},
                        text: [1,2,3]
                        type: 'scatter3d'};

                    var data = [trace2];
                    var layout = {margin: {
                        l: 0,
                        r: 0,
                        b: 0,
                        t: 0
                    }};
                    Plotly.newPlot('myDiv', data, layout);
                    """),callback=false)