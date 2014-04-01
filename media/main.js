
/* Various callbacks */
var updatePHField = function(event, ui) {
    $("#phField").val(ui.value);
};
var updatePMgField = function(event, ui) {
    $("#pmgField").val(ui.value);
};
var updateISField = function(event, ui) {
    $("#ionStrengthField").val(ui.value);
};
var updateERedField = function(event, ui) {
    $("#electronReductionPotentialField").val(ui.value);
};
var toggleCustomConcentrations = function() {
    var conditions = $("form input=radion[name=conditions]:checked");
    var showConcentrations = conditions.val() == "custom";
    $(".customConcentrations").toggle(showConcentrations);
}
    
$(document).ready(function(){
    // Set up autocomplete
    var options = {
      serviceUrl: '/suggest',
      delimiter: /(^|\+|=|<=>|=>)+\s+\d*\s*/
    };
    var queryField = $('#queryField');
    if (queryField) {
        queryField.autocomplete(options);
    }
    
    // Advanced settings sliders.
    var phSlider = $("#phSlider");
    var pmgSlider = $("#pmgSlider");
    var ionStrengthSlider = $("#ionStrengthSlider");
    var electronReductionPotentialSlider = $("#electronReductionPotentialSlider");
    
    if (phSlider) {
        phSlider.slider({
            min: 0.0,
            max: 14.0,
            step: 0.25,
            value: $("#phField").val(),
            slide: updatePHField,
            change: updatePHField});
    }
    if (pmgSlider) {
        pmgSlider.slider({
            min: 0.0,
            max: 14.0,
            step: 0.25,
            value: $("#pmgField").val(),
            slide: updatePMgField,
            change: updatePMgField});
    }
    if (ionStrengthSlider) {
        ionStrengthSlider.slider({
            min: 0.0,
            max: 0.5,
            step: 0.025,
            value: $("#ionStrengthField").val(),
            slide: updateISField,
            change: updateISField});    
    }
    if (electronReductionPotentialSlider) {
        electronReductionPotentialSlider.slider({
            min: -2000,
            max: 2000,
            step: 20,
            value: $("#electronReductionPotentialField").val(),
            slide: updateERedField,
            change: updateERedField});    
    }
    
    var customConcRadio = $("#customConcRadio");
    if (customConcRadio) {
        toggleCustomConcentrations();
    }
    
    var rxnForm =  $("#rxnForm");
    if (rxnForm) {
        rxnForm.change(toggleCustomConcentrations);        
    }

    // Enable lightbox where desired.
    $('a.lightBox').lightBox();

    // Enable buttons where desired.
    $('.buttonSet').buttonset();
});
