<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xmlns:py="http://purl.org/kid/ns#"
    py:extends="'master.kid'">
<head>
<script src="/Deliverator/static/javascript/myForm.js"/>
<script>
function SelectAllList(CONTROL){
    for (var i = 0; i &lt; CONTROL.length; i++) {
        CONTROL.options[i].selected = true;
    }
}
function DeSelectAllList(CONTROL){
    for (var i = 0; i &lt; CONTROL.length; i++) {
        CONTROL.options[i].selected = false;
    }
}
function update_count(){
    var url = '/Deliverator/countImages';
    var params = {
        Width: myFormContents('form_Width'),
        Unit: myFormContents('form_Unit'),
        Field1: myFormContents('form_Field1'),
        Field2: myFormContents('form_Field2'),
        Field3: myFormContents('form_Field3'),
        parameterfileID: myFormContents('form_ParameterFile'),
        Axis: myFormContents('form_Axis'),
        Type: myFormContents('form_Type'),
        enzorunID: myFormContents('form_enzorunID')
    };
    var data = queryString(params);

    function update_count_callback(result){
        result = evalJSONRequest(result);
        getElement('number_of_images').innerHTML = result;
    }

var req = getXMLHttpRequest();
req.open('POST', url, true);
req.setRequestHeader('Content-Type', 'application/x-www-form-urlencoded;charset=utf-8;'); 
sendXMLHttpRequest(req,data).addCallback(update_count_callback);

}
</script>
<meta content="text/html; charset=utf-8" http-equiv="Content-Type" py:replace="''"/>
<title>The Deliverator - Select Images</title>
</head>
<body onload="js:update_count();">
  <div id="status_block">Select your criteria</div>
  <div id="sidebar">
    <p><span id="number_of_images">0</span> images will be returned</p>
    <p><a href="listParams?rid=${rid}">List Parameter Files</a>
    </p>
  </div>
  <div id="getting_started">
    <?python
    attrs={'rid':rid}
    ?>
    <div id="formselect">
      ${form(action=action, attrs=attrs)}
    </div>
  </div>
</body>
</html>
