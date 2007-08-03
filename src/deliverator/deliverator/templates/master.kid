<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<?python import sitetemplate ?>
<html xmlns="http://www.w3.org/1999/xhtml" xmlns:py="http://purl.org/kid/ns#" py:extends="sitetemplate">

<head py:match="item.tag=='{http://www.w3.org/1999/xhtml}head'" py:attrs="item.items()">
    <meta content="text/html; charset=UTF-8" http-equiv="content-type" py:replace="''"/>
    <title py:replace="''">Your title goes here</title>
    <meta py:replace="item[:]"/>
    <style type="text/css">
        #pageLogin
        {
            font-size: 10px;
            font-family: verdana;
            text-align: right;
        }
    </style>
    <style type="text/css" media="screen">
@import "/Deliverator/static/css/style.css";
</style>
<script src="/Deliverator${tg.tg_js}/MochiKit.js"/>
</head>

<body py:match="item.tag=='{http://www.w3.org/1999/xhtml}body'" py:attrs="item.items()">
    <div id="tophalf">
      <div id="sidebarpadder"> </div>
      <div id="header"> </div>
    </div>
    <div id="sitesidebar">
      <h2>Stuff to do</h2>
      <ul class="links">
        <li><a href="/Deliverator/">Home</a></li>
        <li><a href="selectRun">Select a Run</a></li>
        <li py:if="tg.identity.anonymous">
          <a href="/Deliverator/login">Login</a>
        </li>
        <li py:if="not tg.identity.anonymous">
          <a href="/Deliverator/logout">Logout</a>
        </li>
      </ul>
      <div id="about">
        <a href="http://turbogears.org/" border="0">
          <img border="0" height="48" align="middle" 
            src="/Deliverator/static/images/tg_under_the_hood.png"
            alt="TurboGears under the hood" />
        </a>
        <p><a href="http://yt.spacepope.org/">yt</a> and the Deliverator by Matt Turk</p>
      </div>
    </div>
    <div id="main_content">
    <div py:if="tg_flash" class="flash" py:content="tg_flash"></div>

    <div py:replace="[item.text]+item[:]"/>

	<!-- End of main_content -->
	</div>
</body>

</html>
