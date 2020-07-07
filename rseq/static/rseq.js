$(document).ready(function() {
  var table = null;
  $.getJSON('/_get_table', function(data) {
    $("#run_table").html(data.run_table);
    table = $("#run_table").DataTable();
  });
});

var data = [
  ["", "Ford", "Tesla", "Toyota", "Honda"],
  ["2017", 10, 11, 12, 13],
  ["2018", 20, 11, 14, 13],
  ["2019", 30, 15, 12, 13]
];

var container = document.getElementById('example');
var hot = new Handsontable(container, {
  data: data,
  rowHeaders: true,
  colHeaders: true,
  filters: true,
  dropdownMenu: true
});

