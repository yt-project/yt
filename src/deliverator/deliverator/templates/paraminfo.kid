<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xmlns:py="http://purl.org/kid/ns#"
    py:extends="'master.kid'">
<head>
<?python
import time, types, cPickle
?>
<meta content="text/html; charset=utf-8" http-equiv="Content-Type" py:replace="''"/>
<title>Image Gallery</title>
</head>
<body>
  <div class="gallery_header">
    Selection criteria returned ${len(pfs)}
  </div>
  <table class="paraminfo">
    <tr py:for="p in pfs">
<?python
st = time.ctime(p.GeneratedAt)
pc = cPickle.loads(str(p.EnzoHierarchy))
?>
    <td>
    <table align="top">
    <tr><td valign="top"><b>Generated At:</b></td><td valign="top">${st}</td></tr>
    <tr><td valign="top"><b>Filename:</b></td><td valign="top">${p.FileName}</td></tr>
    <tr><td valign="top"><b>MetaData:</b></td><td valign="top"><pre>${p.metaData}</pre></td></tr>
    <tr py:if="p.enzorun.user == myuser"><td valign="bottom"><a href="deleteParamFiles?toDelete=${p.GeneratedAt}">Delete All Images</a></td></tr>
    </table>
    </td>
    
    <td py:if="isinstance(pc, types.DictType)">
    <table>
    <tr py:for="k in pc.keys()">
    <td valign="top"><b>${k}:</b></td>
    <td valign="top">${pc[k]} </td>
    </tr>
    </table>
    </td>
    </tr>
  </table>
</body>
</html>
