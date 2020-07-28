$(document).ready(function() {
  var table = null;
  $.getJSON('/_get_table', function(data) {
    $("#run_table").html(data.run_table);
    table = $("#run_table").DataTable({
      "order": [[ 0, "desc" ]]
    });
  });
});

// from https://stackoverflow.com/questions/46995456/showing-a-bootstrap-alert-message-for-a-specific-time-in-html-and-css
$(function () {
    var duration = 4000; // 4 seconds
    setTimeout(function () { $('#mainAlertMessage').fadeOut(); }, duration);
});

// from https://stackoverflow.com/questions/48613992/bootstrap-4-file-input-doesnt-show-the-file-name
$('#sample_sheet').on('change',function(){
    //get the file name
    var fileName = $(this).val().replace('C:\\fakepath\\', " ");
    //replace the "Choose a file" label
    $(this).next('.custom-file-label').html(fileName);
})

