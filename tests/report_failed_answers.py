"""
Script to generate report of failed answer tests or to generate golden answers
on cloud platforms like Travis

"""


import argparse
import base64
import collections
import datetime
import logging
import os
import re
import shutil
import tempfile
import xml.etree.ElementTree as ET

import nose
import numpy
import requests

from yt.config import ytcfg
from yt.utilities.answer_testing.framework import AnswerTesting
from yt.utilities.command_line import FileStreamer

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("yt_report_failed_answers")
numpy.set_printoptions(threshold=5, edgeitems=1, precision=4)


def generate_failed_answers_html(failed_answers):
    """Generates html for the failed answer tests

    This function creates a html and embeds the images (actual, expected,
    difference) in it for the failed answers.

    Parameters
    ----------
    failed_answers : dict mapping string to dict
        the key is a string denoting the test name, the value is a
        dictionary that stores the actual, expected and difference plot
        file locations of the test.

    Returns
    -------
    string
        a html page

    """

    html_template = """
    <html><head>
    <style media="screen" type="text/css">
    img{{
      width:100%;
      max-width:800px;
    }}
    </style>
    <h1 style="text-align: center;">Failed Answer Tests</h1>
    <p>
      This report shows images of answer tests that failed when running
      the answer tests.
    </p>
    <p>
      <strong>Acutal Image:</strong> plot generated while running the test<br/>
      <strong>Expected Image:</strong> golden answer image<br/>
      <strong>Difference Image:</strong> difference in the "actual"
      and "expected" image
    </p>
    <hr/>
    </head><body>
    <table>{rows}</table>
    </body></html>
    """

    row_template = """
    <tr>
    <td align="center">Actual</td>
    <td align="center">Expected</td>
    <td align="center">Difference</td>
    </tr>
    <tr>
    <td><img src="data:image/png;base64,{0}"></td>
    <td><img src="data:image/png;base64,{1}"></td>
    <td><img src="data:image/png;base64,{2}"></td>
    </tr>
    <tr><td align="center" colspan="3"><b>Test: {3}</b><hr/></td></tr>
    """

    rows = []

    for failed_test_file in failed_answers.values():
        for test_name, images in failed_test_file.items():
            encoded_images = {}
            for key in images:
                with open(images[key], "rb") as img:
                    img_data = base64.b64encode(img.read()).decode()
                    encoded_images[key] = img_data

            formatted_row = row_template.format(
                encoded_images["Actual"],
                encoded_images["Expected"],
                encoded_images["Difference"],
                test_name,
            )
            rows.append(formatted_row)

    html = html_template.format(rows="\n".join(rows))
    return html


def upload_to_curldrop(data, filename):
    """Uploads file to yt's curldrop server

    Uploads bytes `data` by the name `filename` to yt curldrop server.

    Parameters
    ----------
    data : bytes
        Content to be uploaded

    filename : string
        Name of file at curldrop's upload server

    Returns
    -------
    requests.models.Response
        Response returned by curldrop server

    """
    if "TRAVIS" in os.environ:
        job_num = os.environ["TRAVIS_JOB_NUMBER"]
        file_id = "Travis_Job_Num_" + job_num.replace(".", "_")
    else:
        file_id = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = filename.format(file_id)

    base_url = ytcfg.get("yt", "curldrop_upload_url")
    upload_url = base_url + "/" + os.path.basename(filename)
    response = requests.put(upload_url, data=data)
    return response


def upload_failed_answers(failed_answers):
    """Uploads the result of failed answer tests

    Uploads a html page of the failed answer tests.

    Parameters
    ----------
    failed_answers : dict mapping string to dict
        the key is a string denoting the test name, the value is a
        dictionary that stores the actual, expected and difference plot
        file locations of the test.

    Returns
    -------
    requests.models.Response
        Response as returned by the upload service

    """
    html = generate_failed_answers_html(failed_answers)
    # convert html str to bytes
    html = html.encode()
    response = upload_to_curldrop(data=html, filename="failed_answers_{}.html")

    return response


def generate_answers(answer_dir, answers):
    """Generate golden answers

    Generates golden answers for the list of answers in ``answers`` and
    saves them at ``answer_dir``.

    Parameters
    ----------
    answer_dir : string
        directory location to save the generated answers

    answers : list of string
        Collection of missing answer tests specifying full name of the test.
        eg. ['yt.visualization.tests.test_line_plots:test_multi_line_plot']

    Returns
    -------
    bool
        True, if all the missing answers are successfully generated
        False, otherwise

    """
    status = True
    test_argv = [
        os.path.basename(__file__),
        "--with-answer-testing",
        "--nologcapture",
        "-s",
        "-d",
        "-v",
        "--local",
        f"--local-dir={answer_dir}",
        "--answer-store",
    ]

    for job in answers:
        log.info("\n Generating answers for %s", job)
        status &= nose.run(
            argv=test_argv + [job], addplugins=[AnswerTesting()], exit=False
        )
    return status


def upload_answers(answers):
    """Uploads answers not present in answer-store

    This function generates the answers for tests that are not present in
    answer store and uploads a zip file of the same.

    Parameters
    ----------
    answers : list of string
        Collection of missing answer tests specifying full name of the test.
        eg. ['yt.visualization.tests.test_line_plots:test_multi_line_plot']

    Returns
    -------
    requests.models.Response
        Response as returned by the upload service when answers are
        successfully uploaded

    None
        for the case when there was some error while generating the missing
        golden-answers

    """
    # Create temporary location to save new answers
    tmpdir = tempfile.mkdtemp()
    answer_dir = os.path.join(tmpdir, "answer-store")
    if not os.path.exists(answer_dir):
        os.mkdir(answer_dir)
    zip_file = os.path.join(tmpdir, "new-answers")

    status = generate_answers(answer_dir, answers)
    if status:
        zip_file = shutil.make_archive(zip_file, "zip", answer_dir)
        data = iter(FileStreamer(open(zip_file, "rb")))
        response = upload_to_curldrop(data=data, filename="new_answers_{}.zip")
        shutil.rmtree(tmpdir)
        return response
    return None


def extract_image_locations(error_string):
    """Regex based function to extract image file locations.

    Parameters
    ----------
    error_string : String
        The input string having file locations of 'Actual', 'Expected' and
        'Difference' plots. This string is generated by yt's answer-testing
        plugin, when the plot generated in the test does not match to its
        golden answer image.

    Returns
    -------
    dict
        If the `error_string` is successfully parsed to extract plot locations,
        then a dictionary with the keys 'Actual', 'Expected','Difference' and
        values having corresponding plot file locations is returned.
        eg. {'Actual': '/usr/tmp/tmp43la9b0w.png',
             'Expected': '/usr/tmp/tmpbpaqbgi3.png',
             'Difference': '/usr/tmp/tmp43la9b0w-failed-diff.png'}
    None
        When `error_string` does not conform to yt's answer-testing error
        message, which has the information for plot file locations on disk.

    """
    unknown_failure = False
    base_regex = r"\s*\n\s*(.*?.png)"
    img_regex = {
        "Actual": "Actual:" + base_regex,
        "Expected": "Expected:" + base_regex,
        "Difference": "Difference:" + base_regex,
    }
    img_path = {}
    for key in img_regex:
        result = re.search(img_regex[key], error_string, re.MULTILINE)
        if not result:
            unknown_failure = True
            break
        # store the locations of actual, expected and diff plot files
        img_path[key] = result.group(1)

    if not unknown_failure:
        return img_path
    return None


def parse_nose_xml(nose_xml):
    """Parse xml file generated by nosetests.

    Parse nose xml file to find following details:
        Failed tests: These could be due to difference in golden answer image
        and corresponding test plot.

        Missing tests: These errors occur when a corresponding golden answer
        image is not found.

    Parameters
    ----------
    nose_xml : string
        full path of xml file to be parsed

    Returns
    -------
    tuple : (failed_answers, missing_answers)

        failed_answers : list of tuples (string, dict)
        Collection of tuples where the first part is a string denoting the
        test name, the second part is a dictionary that stores the actual,
        expected and difference plot file locations of the test.
        eg. [('yt.visualization.tests.test_line_plots:test_line_plot',
                {'Actual': '/usr/tmp/tmp43la9b0w.png',
                'Expected': '/usr/tmp/tmpbpaqbgi3.png',
                'Difference': '/usr/tmp/tmp43la9b0w-failed-diff.png'}
            )]

        missing_answers : list of string
        Collection of missing answer tests specifying full name of the test.
        eg. ['yt.visualization.tests.test_line_plots:test_multi_line_plot']

    """
    missing_answers = set()
    failed_answers = collections.defaultdict(lambda: dict())
    missing_errors = ["No such file or directory", "There is no old answer available"]
    tree = ET.parse(nose_xml)
    testsuite = tree.getroot()

    for testcase in testsuite:
        for error in testcase.iter("error"):
            handle_error(
                error, testcase, missing_errors, missing_answers, failed_answers
            )
        for error in testcase.iter("failure"):
            handle_error(
                error, testcase, missing_errors, missing_answers, failed_answers
            )
    return failed_answers, missing_answers


def handle_error(error, testcase, missing_errors, missing_answers, failed_answers):
    attribs = ["classname", "name"]
    test_name = ":".join(testcase.attrib[a] for a in attribs)
    message = error.attrib["message"]
    if (
        missing_errors[0] in error.attrib["message"]
        or missing_errors[1] in error.attrib["message"]
    ):
        missing_answers.add(test_name)
    elif "Items are not equal" in error.attrib["message"]:
        img_path = extract_image_locations(error.attrib["message"])
        if img_path:
            failed_answers[test_name][message] = img_path


if __name__ == "__main__":
    """Report failed answer tests of cloud platforms like Travis, Appveyor

    This script parses the nosetests xml file generated after answer tests are
    executed. If the test fail due to difference in actual and expected images,
    this function uploads a html page having all the plots which got failed
    (if executed with `-f` command line argument).
    In case, answer store does not has a golden answer and if executed with
    `-m` argument, it uploads missing answers zip file to yt's curldrop server.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--upload-failed-tests",
        action="store_true",
        help="Upload a comparison report of failed answer tests"
        " to yt's curldrop server.",
    )
    parser.add_argument(
        "-m",
        "--upload-missing-answers",
        action="store_true",
        help="Upload tests' answers that are not found in answer-store.",
    )
    parser.add_argument(
        "--xunit-file",
        action="store",
        dest="nosetest_xml",
        required=True,
        help="Name of the nosetests xml file to parse for failed answer tests.",
    )
    args = parser.parse_args()

    # ANSI color codes
    COLOR_PURPLE = "\x1b[35;1m"
    COLOR_CYAN = "\x1b[36;1m"
    COLOR_RESET = "\x1b[0m"
    FLAG_EMOJI = " \U0001F6A9 "

    failed_answers = missing_answers = None
    if args.upload_failed_tests or args.upload_missing_answers:
        failed_answers, missing_answers = parse_nose_xml(args.nosetest_xml)

    if args.upload_failed_tests and failed_answers:
        response = upload_failed_answers(failed_answers)
        msg = ""
        if response.ok:
            msg += (
                "\n"
                + FLAG_EMOJI
                + COLOR_PURPLE
                + "Successfully uploaded failed answer test(s) result."
                " More details about the test failure can be found at the"
                " URL: "
                + response.text.split("\n")[1]
                + COLOR_RESET
                + FLAG_EMOJI
                + "\n"
            )
        response = upload_answers(failed_answers)
        if response.ok:
            msg += (
                FLAG_EMOJI
                + COLOR_CYAN
                + "Successfully uploaded answer(s) for failed test at URL: "
                + response.text.split("\n")[1]
                + " . Please commit these "
                "answers in the repository's answer-store." + COLOR_RESET + FLAG_EMOJI
            )
            log.info(msg)

    if args.upload_missing_answers and missing_answers:
        response = upload_answers(missing_answers)
        if response.ok:
            msg = (
                FLAG_EMOJI
                + COLOR_CYAN
                + "Successfully uploaded missing answer(s) at URL: "
                + response.text.split("\n")[1]
                + " . Please commit these "
                "answers in the repository's answer-store." + COLOR_RESET + FLAG_EMOJI
            )
            log.info(msg)
