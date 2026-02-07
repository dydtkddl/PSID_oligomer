#!perl

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

my %Args;
GetOptions(\%Args, "Structure=s", "JobDir=s", "Forcefield=s", "DoRefine=s");
my $xtd = Documents->Import("$Args{JobDir}/$Args{Structure}.arc");
my $trj = $xtd->Trajectory;
my $xtdNew = Documents->New($xtd->Name . "_.xtd");
my $trjNew = $xtdNew->Trajectory;

my $currentForcefield = $Args{Forcefield};
$currentForcefield = "COMPASSIII" if($currentForcefield eq "compass");

Modules->Forcite->ChangeSettings([
   CurrentForcefield => $currentForcefield,
   "3DPeriodicvdWSummationMethod" => "Atom based",
   "3DPeriodicvdWChargeGroupCubicSplineCutOff" => 8.5,
   "3DPeriodicvdWChargeGroupCubicSplineWidth" => 0.0,
   "3DPeriodicNonBond2BodyListBufferWidth" => 0.5,
   "3DPeriodicElectrostaticSummationMethod" => "Ewald",
   WriteLevel => "Silent"
]);

my $doRefine = $Args{DoRefine} =~ /[YyTt]/;
my $dynamics = Modules->Forcite->Dynamics;
my $dynamicsSettings = Settings(
   NumberOfSteps => 1000,
   Temperature => 298,
   Ensemble3D => "NVT",
   Thermostat => "Velocity Scale",
);
my $failed;
my $numFrames = $trj->NumFrames;
for(my $frame = 1; $frame <= $numFrames && !defined $failed; $frame++)
{
   $trj->CurrentFrame = $frame;
   my $doc = Documents->New("temp.xsd");
   $doc->CopyFrom($xtd);
   if ( $doRefine )
   {
       eval
       {
          print "Refining frame $frame of $numFrames...";
          $dynamics->Run($doc,$dynamicsSettings)->Trajectory->Delete;
       };
       $failed = ($@) ? $@ : undef;
       my $status = (defined $failed) ? "failed" : "OK";
       print "$status\n";
   }
   $trjNew->AppendFramesFrom($doc);
   $doc->Discard;
}
$xtdNew->Export("$Args{JobDir}/$Args{Structure}.xtd");
die $failed if defined $failed;
print "Refinement finished successfully" if $doRefine;
