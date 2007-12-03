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
<meta content="text/html; charset=utf-8" http-equiv="Content-Type" py:replace="''"/>
<title>The Deliverator - Image Gallery</title>
</head>
<body>
  <div id="status_block">Selection criteria returned ${len(images)}</div>
<?python
NUM = 3
k = []
for i, image in enumerate(images):
    if i % NUM == 0:
        k.append([])
    k[-1].append(image)
splitImages = k
?>
  <form action="deleteImages" method="post">
  <table class="imagegallery">
    <tr py:for="row in splitImages">
    <td py:for="image in row" class="imagebox">
        <a href="viewimage?id=${image.id}">
        <img src="${image.IMG_src}" width="250"/>
        </a>
        <div id="imagecaption">
        <span py:if="show_delete == True">
            <INPUT TYPE="CHECKBOX" NAME="todelete" VALUE="${image.id}">Delete Entry</INPUT>
        </span>
        <ul>
        <li><b>ParameterFile:</b>
            <a href="paraminfo?id=${image.parameterfile.GeneratedAt}">
            ${image.parameterfile.FileName.split("/")[-1]}</a></li>
        <li><b>Type:</b> ${image.Type}</li>
        <li><b>Axis:</b> ${image.Axis}</li>
        <li><b>Width:</b> ${image.Width}</li>
        <li><b>Unit:</b> ${image.Unit}</li>
        <li><b>Field1:</b> ${image.Field1}</li>
        <li><b>Field2:</b> ${image.Field2}</li>
        <li><b>Field3:</b> ${image.Field3}</li>
        <li>(<a href="${image.IMG_src}">view raw</a>)</li>
        </ul>
        </div>
    </td>
    </tr>
  </table>
  <INPUT py:if="show_delete == True" TYPE="SUBMIT" NAME="Delete" VALUE="Delete"/>
  <span py:if="show_delete == True"><input type="checkbox" value="on" name="allbox" onclick="checkAll();"/>Check All</span>
  </form>
</body>
</html>
