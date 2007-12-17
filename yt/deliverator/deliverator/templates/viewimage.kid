<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xmlns:py="http://purl.org/kid/ns#"
    py:extends="'master.kid'">
<head>
<meta content="text/html; charset=utf-8" http-equiv="Content-Type" py:replace="''"/>
<title>The Deliverator - Image Gallery</title>
</head>
<body>
  <div id="status_block">Image View (${Image.id})</div>
  <div id="imageview">
    <a href="${Image.IMG_src}">
    <img src="${Image.IMG_src}" width="500"/>
    </a>
    <div id="Imagecaption">
    <ul>
    <li><b>ParameterFile:</b>
      <a href="paraminfo?id=${Image.parameterfile.GeneratedAt}">
      ${Image.parameterfile.FileName.split("/")[-1]}</a></li>
    <li><b>Type:</b> ${Image.Type}</li>
    <li><b>Axis:</b> ${Image.Axis}</li>
    <li><b>Width:</b> ${Image.Width}</li>
    <li><b>Unit:</b> ${Image.Unit}</li>
    <li><b>Field1:</b> ${Image.Field1}</li>
    <li><b>Field2:</b> ${Image.Field2}</li>
    <li><b>Field3:</b> ${Image.Field3}</li>
    </ul>
    </div>
  </div>
</body>
</html>
