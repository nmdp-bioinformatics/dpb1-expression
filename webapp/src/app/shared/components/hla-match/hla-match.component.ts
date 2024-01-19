import { Component, OnInit, Input } from '@angular/core';
import { Subject } from '@app/shared/models/subject/subject.model';

@Component({
  selector: 'app-hla-match',
  templateUrl: './hla-match.component.html',
  styleUrls: ['./hla-match.component.scss']
})
export class HlaMatchComponent implements OnInit {
  @Input() subject: Subject;

  constructor() { }

  ngOnInit() {
  }

}
