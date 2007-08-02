<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xmlns:py="http://purl.org/kid/ns#"
    py:extends="'master.kid'">
<head>
<meta content="text/html; charset=utf-8" http-equiv="Content-Type" py:replace="''"/>
<title>Welcome to TurboGears</title>
</head>
<body>

  <div id="status_block">The Deliverator is now running</div>
  <div id="sidebar">
    <h2>Places to go</h2>
    <ul class="links">
     <li><a href="login">Login</a></li>
     <li><a href="selectRun">Select a Run</a></li>
    </ul>
  </div>
  <div id="getting_started">
    <ol id="getting_started_steps">
      <li class="getting_started">
        <h3>YT</h3>
        <p>Use <span class="code"><a href="http://www.stanford.edu/~mturk/yt.html">yt</a></span><br/>
        If you're using <span class="code">red</span>, there's an install in the directory <span class="code">/usr/work/mturk/local/</span>.
        </p>
      </li>
      <li class="getting_started">
        <h3>Get an API Key and username</h3>
        <p>Email <a href="mailto:mturk@slac.stanford.edu">Matt</a> for an API
            key and Deliverator username.  Place the API key in
            <span class="code">~/.yt/config</span> in the format:
            <br/><br/><span class="code">[Deliverator]<br/>api-key: YOURKEY</span><br/><br/>
            and use your login on this website.
        </p>
      </li>
      <li class="getting_started">
        <h3>Submit!</h3>
        <p>Use the <span class="code">yt.deliverator.upload</span> module to
            upload image locations to The Deliverator.  See examples in the yt
            distribution for more information.  (Check out the <span class="code">test_fido.py</span> script.)
        </p>
      </li>
    </ol>
  </div>
  <!-- End of getting_started -->
</body>
</html>
