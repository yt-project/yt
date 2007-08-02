<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xmlns:py="http://purl.org/kid/ns#"
    py:extends="'master.kid'">
<head>
<script>
function checkAll(){
    for (var i=0;i&lt;document.forms[0].elements.length;i++)
    {
        var e=document.forms[0].elements[i];
        if ((e.name != 'allbox') &amp;&amp; (e.type=='checkbox'))
        {
            e.checked=document.forms[0].allbox.checked;
        }
    }
}
</script>
<?python
import time,os.path
?>
<meta content="text/html; charset=utf-8" http-equiv="Content-Type" py:replace="''"/>
<title>Parameter Files</title>
</head>
<body>
  <div class="gallery_header">
    Selection criteria returned ${len(run.ParameterFiles)}
  </div>
    <form action="deleteParamFiles" method="post">
    <table align="top">
    <tr><th>Generated At</th><th>FileName</th></tr>
    <tr py:for="p in run.ParameterFiles">
    <td valign="top">${time.ctime(p.GeneratedAt)}</td>
    <td valign="top"><a href="paraminfo?id=${p.GeneratedAt}">${os.path.basename(p.FileName)}</a></td>
    <td py:if="show_delete == True">
        <INPUT TYPE="CHECKBOX" NAME="todelete" VALUE="${p.GeneratedAt}">Delete Entry</INPUT>
    </td>
    </tr>
    </table>
    <INPUT py:if="show_delete == True" TYPE="SUBMIT" NAME="Delete" VALUE="Delete"/>
    <span py:if="show_delete == True"><input type="checkbox" value="on" name="allbox" onclick="checkAll();"/>Check All</span>
    </form>
    <div py:if="show_delete==True" id="status_block"><b><a href="deleteRun?todelete=${run.id}">Delete Run</a></b>
    </div>
</body>
</html>
