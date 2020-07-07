$(document).ready(function() {
  var table = null;
  $.getJSON('/_get_table', function(data) {
    $("#run_table").html(data.run_table);
    table = $("#run_table").DataTable();
  });
});




