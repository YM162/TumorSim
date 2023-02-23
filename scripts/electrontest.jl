using DrWatson
@quickactivate "TumorSim"

using TumorSim
using Blink

w = Window()
body!(w,String(read("scripts/probando.html")),async=false)
loadjs!(w,"https://cdn.plot.ly/plotly-2.18.0.min.js")

js(w, Blink.JSString("""TESTER = document.getElementById('tester')"""))

js(w, Blink.JSString("""Plotly.newPlot( TESTER, [{x: [1, 2, 3, 4, 5],y: [1, 2, 4, 8, 16] }], {margin: { t: 0 } } );"""),callback=false)