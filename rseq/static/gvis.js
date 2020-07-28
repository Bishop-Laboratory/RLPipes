var d3 = require("d3");


const myScript = document.getElementById('bundle_gvis');

var digraph = myScript.attr('digraph')
var elemname = myScript.attr('elemname')

console.log(digraph)

d3.select("#canvas")
  .graphviz()
    .renderDot('digraph {a -> b}');

