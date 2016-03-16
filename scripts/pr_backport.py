from __future__ import print_function
import hglib
import requests
import shutil
import tempfile

from datetime import datetime
from distutils.version import LooseVersion
from time import strptime, mktime

MERGED_PR_ENDPOINT = ("http://bitbucket.org/api/2.0/repositories/yt_analysis/"
                      "yt/pullrequests/?state=MERGED")

YT_REPO = "https://bitbucket.org/yt_analysis/yt"


def clone_new_repo(source=None):
    """Clones a new copy of yt_analysis/yt and returns a path to it"""
    path = tempfile.mkdtemp()
    dest_repo_path = path+'/yt-backport'
    if source is None:
        source = YT_REPO
    hglib.clone(source=source, dest=dest_repo_path)
    with hglib.open(dest_repo_path) as client:
        # Changesets that are on the yt branch but aren't topological ancestors
        # of whichever changeset the experimental bookmark is pointing at
        bookmarks, _ = client.bookmarks()
        bookmark_names = [b[0] for b in bookmarks]
        if 'experimental' in bookmark_names:
            client.update('heads(branch(yt) - ::bookmark(experimental))')
        else:
            client.update('heads(branch(yt))')
    return dest_repo_path


def get_first_commit_after_last_major_release(repo_path):
    """Returns the SHA1 hash of the first commit to the yt branch that wasn't
    included in the last tagged release.
    """
    with hglib.open(repo_path) as client:
        tags = client.log("reverse(tag())")
        tags = sorted([LooseVersion(t[2]) for t in tags])
        for t in tags[::-1]:
            if t.version[0:2] != ['yt', '-']:
                continue
            if len(t.version) == 4 or t.version[4] == 0:
                last_major_tag = t
                break
        last_before_release = client.log(
            "last(ancestors(%s) and branch(yt))" % str(last_major_tag))
        first_after_release = client.log(
            "first(descendants(%s) and branch(yt) and not %s)"
            % (last_before_release[0][1], last_before_release[0][1]))
    return str(first_after_release[0][1][:12])


def get_branch_tip(repo_path, branch, exclude=None):
    """Returns the SHA1 hash of the most recent commit on the given branch"""
    revset = "head() and branch(%s)" % branch
    with hglib.open(repo_path) as client:
        if exclude is not None:
            try:
                client.log(exclude)
                revset += "and not %s" % exclude
            except hglib.error.CommandError:
                pass
        change = client.log(revset)[0][1][:12]
    return change


def get_lineage_between_release_and_tip(repo_path, first, last):
    """Returns the lineage of changesets that were at one point the public tip"""
    with hglib.open(repo_path) as client:
        lineage = client.log("'%s'::'%s' and p1('%s'::'%s') + '%s'"
                             % (first, last, first, last, last))
        return lineage


def get_pull_requests_since_last_release(repo_path):
    """Returns a list of pull requests made since the last tagged release"""
    r = requests.get(MERGED_PR_ENDPOINT)
    done = False
    merged_prs = []
    with hglib.open(repo_path) as client:
        last_tag = client.log("reverse(tag())")[0]
    while not done:
        if r.status_code != 200:
            raise RuntimeError
        data = r.json()
        prs = data['values']
        for pr in prs:
            activity = requests.get(pr['links']['activity']['href']).json()
            merge_date = None
            for action in activity['values']:
                if 'update' in action and action['update']['state'] == 'MERGED':
                    merge_date = action['update']['date']
                    merge_date = merge_date.split('.')[0]
                    timestamp = mktime(strptime(merge_date, "%Y-%m-%dT%H:%M:%S"))
                    merge_date = datetime.fromtimestamp(timestamp)
                    break
            if merge_date is None:
                break
            if merge_date < last_tag[6]:
                done = True
                break
            merged_prs.append(pr)
        r = requests.get(data['next'])
    return merged_prs


def cache_commit_data(prs):
    """Avoid repeated calls to bitbucket API to get the list of commits per PR"""
    commit_data = {}
    for pr in prs:
        data = requests.get(pr['links']['commits']['href']).json()
        if data.keys() == [u'error']:
            # this happens when commits have been stripped, e.g.
            # https://bitbucket.org/yt_analysis/yt/pull-requests/1641
            continue
        done = False
        commits = []
        while not done:
            commits.extend(data['values'])
            if 'next' not in data:
                done = True
            else:
                data = requests.get(data['next']).json()
        commit_data[pr['id']] = commits
    return commit_data


def find_commit_in_prs(needle, commit_data, prs):
    """Finds the commit `needle` PR in the commit_data dictionary

    If found, returns the pr the needle commit is in. If the commit was not
    part of the PRs in the dictionary, returns None.
    """
    for pr_id in commit_data:
        commits = commit_data[pr_id]
        for commit in commits:
            if commit['hash'] == needle[1]:
                pr = [pr for pr in prs if pr['id'] == pr_id][0]
                return pr
    return None


def find_merge_commit_in_prs(needle, prs):
    """Find the merge commit `needle` in the list of `prs`

    If found, returns the pr the merge commit comes from. If not found, return
    None
    """
    for pr in prs[::-1]:
        if pr['merge_commit'] is not None:
            if pr['merge_commit']['hash'] == needle[1][:12]:
                return pr
    return None


def create_commits_to_prs_mapping(linege, prs):
    """create a mapping from commits to the pull requests that the commit is
    part of
    """
    commits_to_prs = {}
    # make a copy of this list to avoid side effects from calling this function
    my_prs = list(prs)
    commit_data = cache_commit_data(my_prs)
    for commit in lineage:
        cset_hash = commit[1]
        message = commit[5]
        if message.startswith('Merged in') and '(pull request #' in message:
            pr = find_merge_commit_in_prs(commit, my_prs)
            if pr is None:
                continue
            commits_to_prs[cset_hash] = pr
            # Since we know this PR won't have another commit associated with it,
            # remove from global list to reduce number of network accesses
            my_prs.remove(commits_to_prs[cset_hash])
        else:
            pr = find_commit_in_prs(commit, commit_data, my_prs)
            commits_to_prs[cset_hash] = pr
    return commits_to_prs


def invert_commits_to_prs_mapping(commits_to_prs):
    """invert the mapping from individual commits to pull requests"""
    inv_map = {}
    for k, v in commits_to_prs.iteritems():
        # can't save v itself in inv_map since it's an unhashable dictionary
        if v is not None:
            created_date = v['created_on'].split('.')[0]
            timestamp = mktime(strptime(created_date, "%Y-%m-%dT%H:%M:%S"))
            created_date = datetime.fromtimestamp(timestamp)
            pr_desc = (v['id'], v['title'], created_date,
                       v['links']['html']['href'], v['description'])
        else:
            pr_desc = None
        inv_map[pr_desc] = inv_map.get(pr_desc, [])
        inv_map[pr_desc].append(k)
    return inv_map


def get_last_descendant(repo_path, commit):
    """get the most recent descendant of a commit"""
    with hglib.open(repo_path) as client:
        com = client.log('last(%s::)' % commit)
    return com[0][1][:12]

def screen_already_backported(repo_path, inv_map):
    with hglib.open(repo_path) as client:
        tags = client.log("reverse(tag())")
        major_tags = [t for t in tags if t[2].endswith('.0')]
        most_recent_major_tag_name = major_tags[0][2]
        lineage = client.log(
            "descendants(%s) and branch(stable)" % most_recent_major_tag_name)
        prs_to_screen = []
        for pr in inv_map:
            for commit in lineage:
                if commit[5].startswith('Backporting PR #%s' % pr[0]):
                    prs_to_screen.append(pr)
        for pr in prs_to_screen:
            del inv_map[pr]
        return inv_map

def commit_already_on_stable(repo_path, commit):
    with hglib.open(repo_path) as client:
        commit_info = client.log(commit)[0]
        most_recent_tag_name = client.log("reverse(tag())")[0][2]
        lineage = client.log(
            "descendants(%s) and branch(stable)" % most_recent_tag_name)
        # if there is a stable commit with the same commit message,
        # it's been grafted
        if any([commit_info[5] == c[5] for c in lineage]):
            return True
        return False

def backport_pr_commits(repo_path, inv_map, last_stable, prs):
    """backports pull requests to the stable branch.

    Accepts a dictionary mapping pull requests to a list of commits that
    are in the pull request.
    """
    pr_list = inv_map.keys()
    pr_list = sorted(pr_list, key=lambda x: x[2])
    for pr_desc in pr_list:
        merge_warn = False
        merge_commits = []
        pr = [pr for pr in prs if pr['id'] == pr_desc[0]][0]
        data = requests.get(pr['links']['commits']['href']).json()
        commits = data['values']
        while 'next' in data:
            data = requests.get(data['next']).json()
            commits.extend(data['values'])
        commits = [com['hash'][:12] for com in commits]
        with hglib.open(repo_path) as client:
            for com in commits:
                if client.log('merge() and %s' % com) != []:
                    merge_warn = True
                    merge_commits.append(com)
        if len(commits) > 1:
            revset = " | ".join(commits)
            revset = '"%s"' % revset
            message = "Backporting PR #%s %s" % \
                (pr['id'], pr['links']['html']['href'])
            dest = get_last_descendant(repo_path, last_stable)
            message = \
                "hg rebase -r %s --keep --collapse -m \"%s\" -d %s\n" % \
                (revset, message, dest)
            message += "hg update stable\n\n"
            if merge_warn is True:
                if len(merge_commits) > 1:
                    merge_commits = ", ".join(merge_commits)
                else:
                    merge_commits = merge_commits[0]
                message += \
                    "WARNING, PULL REQUEST CONTAINS MERGE COMMITS, CONSIDER\n" \
                    "BACKPORTING BY HAND TO AVOID BACKPORTING UNWANTED CHANGES\n"
                message += \
                    "Merge commits are %s\n\n" % merge_commits
        else:
            if commit_already_on_stable(repo_path, commits[0]) is True:
                continue
            message = "hg graft %s\n" % commits[0]
        print("PR #%s\nTitle: %s\nCreated on: %s\nLink: %s\n%s" % pr_desc)
        print("To backport, issue the following command(s):\n")
        print(message)
        raw_input('Press any key to continue')


if __name__ == "__main__":
    print("")
    print("Gathering PR information, this may take a minute.")
    print("Don't worry, yt loves you.")
    print("")
    repo_path = clone_new_repo()
    try:
        last_major_release = get_first_commit_after_last_major_release(repo_path)
        last_dev = get_branch_tip(repo_path, 'yt', 'experimental')
        last_stable = get_branch_tip(repo_path, 'stable')
        lineage = get_lineage_between_release_and_tip(
            repo_path, last_major_release, last_dev)
        prs = get_pull_requests_since_last_release(repo_path)
        commits_to_prs = create_commits_to_prs_mapping(lineage, prs)
        inv_map = invert_commits_to_prs_mapping(commits_to_prs)
        # for now, ignore commits that aren't part of a pull request since
        # the last bugfix release. These are mostly commits in pull requests
        # from before the last bugfix release but might include commits that
        # were pushed directly to the repo.
        del inv_map[None]

        inv_map = screen_already_backported(repo_path, inv_map)
        print("In another terminal window, navigate to the following path:")
        print("%s" % repo_path)
        raw_input("Press any key to continue")
        backport_pr_commits(repo_path, inv_map, last_stable, prs)
        raw_input(
            "Now you need to push your backported changes. The temporary\n"
            "repository currently being used will be deleted as soon as you\n"
            "press any key.")
    finally:
        shutil.rmtree(repo_path)
