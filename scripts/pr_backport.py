import shutil
import tempfile

import dateutil.parser
import git
import requests

API_URL = "https://api.github.com/graphql"

YT_REPO = "https://github.com/yt-project/yt"

PR_QUERY = """
{
  repository(owner: "yt-project", name: "yt") {
    pullRequests(first:100 states:MERGED %s) {
      edges {
        node {
          number
          title
          mergedAt
          author {
            login
          }
          body
          url
        }
      }
      pageInfo {
        endCursor
        hasNextPage
        hasPreviousPage
        startCursor
      }
    }
  }
}
"""


def clone_new_repo(source=None):
    """Clones a new copy of yt_analysis/yt and returns a path to it"""
    path = tempfile.mkdtemp()
    dest_repo_path = path + "/yt-backport"
    if source is None:
        source = YT_REPO
    git.Repo.clone_from(source, dest_repo_path)
    return dest_repo_path


def get_date_of_last_tag(repo_path):
    repo = git.Repo(repo_path)
    tags = sorted(repo.tags, key=lambda t: t.commit.committed_date)
    return tags[-1].commit.committed_date


def get_prs_since_last_release(date, key):
    headers = {"Authorization": "token %s" % key}
    resp = requests.post(url=API_URL, json={"query": PR_QUERY % ""}, headers=headers)
    ret = []
    while True:
        jsr = resp.json()
        cursor = jsr["data"]["repository"]["pullRequests"]["pageInfo"]["endCursor"]
        if cursor is None:
            break
        prs = jsr["data"]["repository"]["pullRequests"]["edges"]
        for pr in prs:
            pr_date = dateutil.parser.parse(pr["node"]["mergedAt"]).timestamp()
            if pr_date > date:
                ret.append(pr["node"])
        resp = requests.post(
            url=API_URL,
            json={"query": PR_QUERY % ('after:"%s"' % cursor)},
            headers=headers,
        )
    return ret


def backport_prs(repo_path, prs):
    for pr in prs:
        print("")
        print("PR %s" % pr["number"])
        print(pr["title"])
        print(pr["author"]["login"])
        print(pr["body"])
        print(pr["url"])
        print("%s.diff" % pr["url"])
        input("Press any key to continue")


if __name__ == "__main__":
    key = input(
        "Please enter your github OAuth API key\n"
        "See the github help for instructions on how to "
        "generate a personal access token.\n>>> "
    )
    print("")
    print("Gathering PR information, this may take a minute.")
    print("Don't worry, yt loves you.")
    print("")
    repo_path = clone_new_repo()
    try:
        date = get_date_of_last_tag(repo_path)
        prs = get_prs_since_last_release(date, key)
        print("In another terminal window, navigate to the following path:")
        print("%s" % repo_path)
        input("Press any key to continue")
        backport_prs(repo_path, prs)
        input(
            "Now you need to push your backported changes. The temporary\n"
            "repository currently being used will be deleted as soon as you\n"
            "press any key."
        )
    finally:
        shutil.rmtree(repo_path)
