$(document).ready(function() {
  var table = null;
  $.getJSON('/_get_table', function(data) {
    $("#run_table").html(data.run_table);
    table = $("#run_table").DataTable();
  });
});

// from https://stackoverflow.com/questions/46995456/showing-a-bootstrap-alert-message-for-a-specific-time-in-html-and-css
$(function () {
    var duration = 4000; // 4 seconds
    setTimeout(function () { $('#mainAlertMessage').fadeOut(); }, duration);
});


