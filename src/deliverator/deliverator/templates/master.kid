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
    <div py:if="tg.config('identity.on',False) and not 'logging_in' in locals()"
        id="pageLogin">
        <span py:if="tg.identity.anonymous">
            <a href="/Deliverator/login">Login</a>
        </span>
        <span py:if="not tg.identity.anonymous">
            Welcome ${tg.identity.user.display_name}.
            <a href="/Deliverator/logout">Logout</a>
        </span>
    </div>
    <div id="header"> </div>
    <div id="main_content">
    <div py:if="tg_flash" class="flash" py:content="tg_flash"></div>

    <div py:replace="[item.text]+item[:]"/>

	<!-- End of main_content -->
	</div>
<div id="footer"> <img src="/static/images/under_the_hood_blue.png" alt="TurboGears under the hood" />
  <p>TurboGears is a open source front-to-back web development
    framework written in Python</p>
  <p>Copyright &copy; 2006 Kevin Dangoor</p>
  <p>Raven component by Matt Turk</p>
</div>
</body>

</html>
